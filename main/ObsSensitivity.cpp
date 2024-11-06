#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>

//Cern ROOT
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "Math/Interpolator.h"
//My Libraries
#include "MathTools.h"
#include "PhyTools.h"
#include "GetObservedEvents.h"


# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds


using namespace std;

void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename);

int main(int argc, char **argv)
{
    //Detector setting 
    double detectorMassMT = 10.0; //Mass in megaton units
    double detectorExposureYears  = 10.0; //Detector Exposure time in sec: One Year
    double exposureMtonYears = detectorMassMT * MTon * detectorExposureYears * years2sec; // Exposure in [Mton*years]
    
    // Detector Position
    double detectorRadius = Rearth;            // Radius at South Pole
    double detectorPos[3] = {0.0, 0.0, -1.0 * detectorRadius};

    // Paths for input and output files
    std::string fluxPath = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/NuFlux/SP_AziAveraged_solmin/spl-nu-20-01-000.d";
    std::string premPath = "/home/dehy0499/OscProb/PremTables/prem_44layers.txt";

     //Pile Set Up---------------------------------------------------------------
    double llvpHeight = 600; // Height of LLVP in km
    double llvpRadius = Rcmb + llvpHeight; //km
    double depthMin = Rcmb-2500;
    double depthMax = llvpRadius+500; // Distance from the center of the Earth
    double llvpDensityContrast = 2.0; // 2% density contrats for LLVP
    std::string llvpShape = "cake";


    // Simulation Configuration

    // True binning
    int zenithBins=60; // # Bins in zenith/cos(zenith)
    int azimuthBins =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int energyBins =60; // bins in energy

    // Reco binning
    int zenithBinsReco=40; // # Bins in zenith/cos(zenith)
    int azimuthBinsReco =100; // # Bins in azimuth (optimal bins are 110 or 22)
    int energyBinsReco =40; // bins in energy

    // Energy [GeV]    
    double enuMin = 2.0; 
    double enuMax = 10.0;

    //Reconstructed Directions

    // Zenith angle [degrees]
    double zenithMin = 180-TMath::ASin( depthMax/Rocean )*(180.0/TMath::Pi()) ; // min 90
    double zenithMax = 180-TMath::ASin( (depthMin)/Rocean )*(180.0/TMath::Pi()) ; // max 180

    double cosZenithMin = cos(zenithMax*TMath::Pi()/180.0);// Cos(zenithMin) - Minimum possible is -1    
    double cosZenithMax = cos(zenithMin*TMath::Pi()/180.0);// Cos(zenithMax) - Maximum possible is 0

    // Azimuth angle [degrees]
    double azimuthMin = -55.0;
    double azimuthMax =  55.0 ;

    // Neutrino flavor (0: nue, 1: numu, 2: nutau)
    int nuflv = 1;


     // File and folder structure configuration
    int exposureYearsLabel = static_cast<int>(detectorMassMT * detectorExposureYears);

    std::string BinLabel = "Obs"+std::to_string(zenithBinsReco)+"Zen"+ std::to_string(azimuthBinsReco)+"Az"+ std::to_string(energyBinsReco)+"Enu_"
                        +"Int"+std::to_string(zenithBins)+"Zen"+ std::to_string(azimuthBins)+"Az"+ std::to_string(energyBins)+"Enu";

    std::string SimLabel = "Simulation_"+llvpShape+"_"+std::to_string(exposureYearsLabel)+"Mton"+"_"+ std::to_string(static_cast<int>(enuMin)) +"-"+ 
                                 std::to_string(static_cast<int>(enuMax))+"GeV_" + std::to_string(static_cast<int>(zenithMin)) +"-"+ 
                                 std::to_string(static_cast<int>(zenithMax))+"Zen_"  + std::to_string(static_cast<int>(azimuthMin)) +"-"+ 
                                 std::to_string(static_cast<int>(azimuthMax))+"Az";
    
    std::string pathToResults = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/SimulationResults/PreliminaryResults/";
    
    std::string pathToChiSquares = pathToResults + "ObsChi2/";
    std::string chiSquareFolderPath = pathToChiSquares + "ObsChi2_" + SimLabel + "/";

    std::string chiSquareFileName = chiSquareFolderPath + "Obs_cth_Chi2LLVP" + std::to_string(exposureYearsLabel) + "myrs_" + llvpShape + 
                                          "nu" + std::to_string(nuflv) + "_" + BinLabel + ".csv";

    // Create output directories if they do not exist
    mkdir(eventFolderPath.c_str(), 0777);
    mkdir(chiSquareFolderPath.c_str(), 0777);




    
   //Standart Earth-------------------------------------------------------------
   AsimovObsSimulation Simulation;

   Simulation.ThePremTable = premPath;
   Simulation.TheHondaTable = fluxPath;
   Simulation.SetDetectorXYZ(detectorPos);
   Simulation.PileInModel = false;

   //Simulation of Standard Earth
   Simulation.SetIntervals(zenithMin,zenithMax,azimuthMin,azimuthMax,enuMin,enuMax);
   Simulation.SetTrueBinning(zenithBins,azimuthBins,energyBins);
   Simulation.SetRecoBinning(zenithBinsReco,azimuthBinsReco,energyBinsReco);
   Simulation.SetExposure(exposureMtonYears);
   Simulation.flvf=nuflv;


   std::vector<TH2D*>  ObservedEventsStandard = Simulation.GetObsEvents3Dcth();


  
    
    std::ofstream SenvData(chiSquareFolderPath); 


    for (int i = 1; i<= 10; ++i)
    {

        Simulation.ThePremTable = premPath;
        Simulation.PileInModel = true;
        Simulation.ShapeOfPile = llvpShape;
        Simulation.ThePileHight = i*100;
        Simulation.ThePileAperture = 45;
        Simulation.ThePileDensityContrast = llvpDensityContrast;
        Simulation.ThePileChemicalContrast = 0.0;
        std::vector<TH2D*>  ObservedEventsAlternative = Simulation.GetObsEvents3Dcth();

        double chi2totth = 0;

        for (int n = 0; n < ObservedEventsAlternative.size(); ++n)
        {
            chi2totth += Get2DChi2( ObservedEventsStandard[n] , ObservedEventsAlternative[n]);
        }

    //std::cout << "Density %: " << rpct << " Chi2: " << chi2tot << std::endl;

        SenvData << i << " " << chi2totth << ", " << zenithBinsReco << ", " << azimuthBinsReco << ", " << energyBinsReco << ", " << Simulation.ThePileHight  <<  std::endl; 

        
    }

    SenvData.close();
   
    return 0;
   


}

void ExportToCSV(TH2D* stdhist, TH2D* althist, TH2D* hist, std::string filename)
{
    std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
    std::string IntFolder = "/SimulationResults/PreliminaryResults/IntEvents/";

    //std::string Evtfname  = EvtFolder.c_str()+chi2namecth;

    //std::string path2file  = NuTomoPath+IntFolder+filename;
   
    std::ofstream outfile(filename);

    // Recorre los bins y guarda el contenido


    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);

            double nexp = stdhist->GetBinContent(i, j);
            double nobs = althist->GetBinContent(i, j); 
            double diff = hist->GetBinContent(i, j);

            outfile << x << "," << y << "," <<  nexp <<" , " << nobs << " , " << diff << std::endl;
        }
    }

    outfile.close();
}


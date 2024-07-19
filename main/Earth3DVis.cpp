#include "Earth3DModel.h"

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

#include "TApplication.h"
#include "TRootCanvas.h"
#include "TCanvas.h"
#include "TF1.h"


//OSCPROB

#include "PremModel.h"

//Some import Earth Radius
# define Rcmb 3480.0 // Core-Mantle-Boundary Radius (km)
# define Rocean 6368.0 // Earth radius in (km)
# define Rearth 6371.0 // Earth radius in (km)
# define Ratm 6386.0 // Radius of the Atmosphere (km)

//Some units
# define mN   1.67492749804E-27  // Nucleons mass (~ Neutron mass) in kg
# define MTon  1E9  //Metric MegaTon in kg
# define years2sec 3.154E7 // Years in Seconds

int main(int argc, char **argv)
{

    double PileHeight = 1000; // Height of LLVP in km
    double PileRadius = Rcmb + PileHeight; //km
    double DepthMin = Rcmb;
    double DepthMax = PileRadius; // Distance from the center of the Earth
    double PileDensityPct = 3.0; // 2% density contrats for LLVP

  double zenmin = 180-TMath::ASin( DepthMax/Rocean )*(180.0/TMath::Pi()) ; // min 90
  double zenmax = 180-TMath::ASin( (DepthMin)/Rocean )*(180.0/TMath::Pi()) ; // max 180

  double th = 170;

  double cth = cos(th*TMath::Pi()/180.0);

  double phi = 0*TMath::Pi()/180.0;

  // OSCPROB

  std::string PremFolder = "/home/dehy0499/OscProb/PremTables/";
  std::string PremName = "prem_44layers.txt";
  std::string path2prem = PremFolder+PremName;

  OscProb::PremModel prem(path2prem);

  prem.FillPath(cth); // Fill paths from PREM model

  std::vector<OscProb::NuPath> paths = prem.GetNuPath();




  //Mine

  TApplication app("app", &argc, argv);

  Earth3DModel test;

  test.SetModel(path2prem);

  test.SetDirection(cth, phi);
  
  test.SetPile( true, "cake", 2.0, 0.0);

  test.SetLayerProp(23,0.0,0.0);
  
  std::vector<std::vector<double>> EarthPath = test.Create3DPath( );


  std::cout << "Path differences" << std::endl;


    for (int i = 0; i < paths.size() ; ++i)
  {
    
    std::cout << i << " " <<paths[i].length - EarthPath[i][0] << " " << paths[i].density - EarthPath[i][1]  << " " << paths[i].zoa - EarthPath[i][2]<< std::endl;

  }


 std::cout << zenmin << " " << zenmax <<  std::endl;

  
  TCanvas * c = test.EarthCanvas;

  //c->Print("./TEST.png");

  c->Modified(); c->Update();

  
  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  app.Run();
  

return 0;

}
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

#include "TGraph.h"
#include "TMultiGraph.h"

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

  std::cout << zenmin << " " << zenmax << std::endl;

  double th = 150;

  //double cth = cos(th*TMath::Pi()/180.0);

  //double cth = -0.8194; //Mid point lower center

  //double cth = -0.95;

  double cth = cos(th*TMath::Pi()/180.0);

  //double cth =  -0.75;

  double phi = 0.0*TMath::Pi()/180.0;

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

  test.PileThickness = 1000;

  test.aWidth = 45;
  test.SetModel(path2prem);

  test.SetDirection(cth, phi);
  
  test.SetPile( true, "cake", 3.0, 0.0);

  test.SetLayerProp(27,0.0,0.0);
  
  std::vector<std::vector<double>> EarthPath = test.Create3DPath( );


  std::cout << "Mine" << std::endl;


  double sumL = 0;

  TMultiGraph * ThePaths = new TMultiGraph();

  TGraph * MyModel = new TGraph(EarthPath.size());

  TGraph * apcModel = new TGraph(paths.size());


 
  std::ofstream ModelPath("ModelPath.csv"); 

  for (int i = 0; i < EarthPath.size() ; ++i)
  {

    sumL += EarthPath[i][0]; 
    std::cout<< i << "  me " << sumL<< " " << std::setprecision(10) <<EarthPath[i][0] << " " <<EarthPath[i][1]  << " " <<EarthPath[i][2]<< std::endl;

    ModelPath << sumL << " " << EarthPath[i][0] << " " << EarthPath[i][1]<< std::endl;

    MyModel->SetPoint(i,sumL,EarthPath[i][1] );

  }

  std::cout << "     " << std::endl;
  std::cout << "OscProb" << std::endl;

  sumL = 0;


  std::ofstream OscProbPath("OscProbPath.csv"); 
  
  for (int i = 0; i < paths.size() ; ++i)
  {

    sumL += paths[i].length ; 
    std::cout << i << " " << sumL << " " << paths[i].length - EarthPath[i][0]  << " " <<  paths[i].density  << " " << paths[i].zoa << std::endl;

    OscProbPath << sumL << " " << paths[i].length << " " << paths[i].density << std::endl;

    apcModel->SetPoint(i,sumL,paths[i].density  );

  }

  OscProbPath.close();
  ModelPath.close();

  ThePaths->Add(MyModel);
  ThePaths->Add(apcModel);

 std::cout << zenmin << " " << zenmax <<  std::endl;

  
  TCanvas * c = test.EarthCanvas;

  //c->Print("./TEST.png");

  c->Modified(); c->Update();

  
  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");


  TCanvas * c1 = new TCanvas();

  ThePaths->Draw("apl");

  c1->Modified(); c1->Update();

  
  TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();
  rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");



  app.Run();
  

return 0;

}
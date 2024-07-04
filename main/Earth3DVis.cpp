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



int main(int argc, char **argv)
{

  double th = 180.0-0;

  double cth = cos(th*TMath::Pi()/180.0);

  double phi = 0;

  // OSCPROB

  std::string PremFile = "prem_44layers.txt";

  std::string model = "/home/dehy0499/OscProb/PremTables/"+PremFile;

  OscProb::PremModel prem(model);

  prem.FillPath(cth); // Fill paths from PREM model

  std::vector<OscProb::NuPath> paths = prem.GetNuPath();




  //Mine

  TApplication app("app", &argc, argv);

  Earth3DModel test;

  test.SetModel(PremFile);

  test.SetDirection(cth, phi);
  
  test.SetPile( false, "pancake" );

  test.SetLayerProp(23,0.0,0.0);
  
  std::vector<std::vector<double>> EarthPath = test.Create3DPath( );


  std::cout << "Paths" << std::endl;


    for (int i = 0; i < paths.size() ; ++i)
  {
    
    std::cout << i << " " <<paths[i].length-EarthPath[i][0]<< " " << paths[i].density-EarthPath[i][1] << " " << paths[i].zoa << std::endl;

  }

  std::cout << " mine  " << std::endl; 
  
    for (int i = 0; i < EarthPath.size() ; ++i)
  {
    
    std::cout << i << " " << EarthPath[i][0]<< " " << EarthPath[i][1] << " " << EarthPath[i][2] << std::endl;

  }

 std::cout << "Note a problem " << std::endl;

  /*
  TCanvas * c = test.EarthCanvas;

  c->Print("./TEST.png");

  c->Modified(); c->Update();

  
  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  app.Run();
  */

return 0;

}
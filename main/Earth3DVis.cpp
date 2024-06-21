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




int main(int argc, char **argv)
{

  TApplication app("app", &argc, argv);

  Earth3DModel test;

  test.SetModel("prem_44layers.txt");

  double th = 140.0;

  double cth = cos(th*TMath::Pi()/180.0);

  double phi = -45;

  test.SetDirection(cth, phi);

  //test.LLVPIdLayers()


  std::vector<int> LLVPSegments {25,26,27,28,29,30,31};

  test.WhichLayersLLVPs = LLVPSegments;

  //test.aWidth = 45;
  

  //void SetLLVPAtt( LLVPSegments, aWidth, 3.0, 0);
  
  test.ActiveHeterogeneity( true );
   
  
  std::vector<std::vector<double>> EarthPath = test.Create3DPath( );

  TCanvas * c = test.EarthCanvas;

  c->Print("./TEST.png");

  c->Modified(); c->Update();

  
  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  app.Run();

  /*
  double TestL = 0  ;

  //double Ltot = -2.0*6371*cos(th);

  for (int i = 0; i < EarthPath.size(); ++i)
  {


    std::cout << EarthPath[i][0] << " " << EarthPath[i][1] << " " << EarthPath[i][2] << " " << std::endl;
    //TestL = TestL + EarthPath[i][0];
    
  }
  */

  /*

  TApplication app("app", &argc, argv);

  TCanvas* c = new TCanvas("c", "Something", 0, 0, 800, 600);
   TF1 *f1 = new TF1("f1","sin(x)", -5, 5);
   f1->SetLineColor(kBlue+1);
   f1->SetTitle("My graph;x; sin(x)");
   f1->Draw();
   c->Modified(); c->Update();

  
  TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
  app.Run();
*/

return 0;

}
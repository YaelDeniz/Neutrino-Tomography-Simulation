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
#include "TH2D.h"


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

void ExportToCSV(TH2D* hist, std::string filename);


int main(int argc, char **argv)
{


  // Simulation Configuration
  int zenithBins = 100;                      // Zenith bins
  int azimuthBins = 100;                     // Azimuth bins


  Earth3DModel Baselines;

  Baselines.PileThickness = 1000;

  Baselines.aWidth = 45;

  std::string llvpShape = "cake";
  double llvpRadius = Rcmb + Baselines.PileThickness;

  double depthMin = Rcmb - 2500;
  double depthMax = llvpRadius + 500;        // Max distance from Earth's center


  // Zenith angle [rad]
  double zenithMin = 180 - TMath::ASin(depthMax / Rearth) * (180.0 / TMath::Pi());
  double zenithMax = 180 - TMath::ASin(depthMin / Rearth) * (180.0 / TMath::Pi());

  double cosZenithMin = cos(zenithMax * TMath::Pi() / 180.0);
  double cosZenithMax = cos(zenithMin * TMath::Pi() / 180.0);

  // Azimuthal angle [rad]
  double azimuthStart = -50;
  double azimuthEnd = 50;
  double dAzimuth = (azimuthEnd - azimuthStart) / (2.0 * azimuthBins);
  double azimuthMin = (azimuthStart - dAzimuth);
  double azimuthMax = (azimuthEnd + dAzimuth);

  TH2D* LLVPBaseLineHIst = new TH2D("LLVPBaseLine", "LLVP Base Line", zenithBins, cosZenithMin, cosZenithMax, azimuthBins, azimuthMin, azimuthMax);

  double cthi = LLVPBaseLineHIst->GetXaxis()->GetBinCenter(1);
  double phij = LLVPBaseLineHIst->GetYaxis()->GetBinCenter(1)*TMath::Pi()/180.0;

  //-----------------------------------------------------------------------------------------------------------------------------------------

    double PileHeight = 1000; // Height of LLVP in km
    double PileRadius = Rcmb + PileHeight; //km
    double DepthMin = Rcmb;
    double DepthMax = PileRadius; // Distance from the center of the Earth
    double PileDensityPct = 30.0; // 2% density contrats for LLVP

  double zenmin = 180-TMath::ASin( DepthMax/Rocean )*(180.0/TMath::Pi()) ; // min 90
  double zenmax = 180-TMath::ASin( (DepthMin)/Rocean )*(180.0/TMath::Pi()) ; // max 180

  std::cout << zenmin << " " << zenmax << std::endl;

  double eta = 35;
  double th = 180-eta;

  //double cth = cos(th*TMath::Pi()/180.0);

  //double cth = -0.8194; //Mid point lower center

  //double cth = -0.95;

  double cth = cos(th*TMath::Pi()/180.0);

  //double cth =  -0.75;

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

  test.PileThickness = 1000;

  test.aWidth = 45;
  test.SetModel(path2prem);

  test.SetDirection(cth, phi);
  
  test.SetPile( true, "cake", 3.0, 0.0);

  test.SetLayerProp(27,0.0,0.0);

  std::vector<std::vector<double>> EarthPath = test.Create3DPath( );

  test.SetPile( true, "cake", 0.0, 0.0);

  std::vector<std::vector<double>> EarthPathPREM = test.Create3DPath( );

  std::cout << "Mine" << std::endl;


  double sumL = 0;
  double sumLup = 0;
  double sumLlow = 0;
  double  sumL_LLVP = 0;


  TMultiGraph * ThePaths = new TMultiGraph();

  TGraph * MyModel = new TGraph(2*EarthPath.size());

  MyModel->SetMarkerStyle(3);

  TGraph * PointsLLVP = new TGraph();

  PointsLLVP->SetMarkerStyle(23);
  PointsLLVP->SetLineStyle(0);

  TGraph * apcModel = new TGraph(2*paths.size());

  apcModel->SetMarkerStyle(20);

  std::vector<double> MyModelSumL;


 
  std::ofstream ModelPath("ModelPath.csv"); 

  for (int i = 0; i < EarthPath.size() ; ++i)
  {

    sumLup += EarthPath[i][0]; 
    sumLlow = sumLup - EarthPath[i][0];   

    std::cout<< i << " " << sumL<< " " << std::setprecision(10) <<EarthPath[i][0] << " " <<EarthPath[i][1]  << " " <<EarthPath[i][2]<< std::endl;

    ModelPath << sumL << " " << EarthPath[i][0] << " " << EarthPath[i][1]<< std::endl;

    MyModelSumL.push_back(sumLlow);
    MyModelSumL.push_back(sumLup);

    //MyModel->SetPoint(i,sumL,EarthPath[i][1] );

    MyModel->SetPoint(2*i,sumLlow,EarthPath[i][1]  );

    MyModel->SetPoint(2*i+1,sumLup,EarthPath[i][1]  );

  }

  std::cout << "     " << std::endl;
  std::cout << "OscProb" << std::endl;

  sumL = 0;

  sumLup = 0; 

  std::ofstream OscProbPath("OscProbPath.csv"); 
  
  for (int i = 0; i < EarthPathPREM.size() ; ++i)
  {


    sumLup += EarthPathPREM[i][0]; 
    sumLlow = sumLup - EarthPathPREM[i][0];   

    std::cout<< i << " " << sumL<< " " << std::setprecision(10) <<EarthPathPREM[i][0] << " " <<EarthPathPREM[i][1]  << " " <<EarthPathPREM[i][2]<< std::endl;


    //MyModel->SetPoint(i,sumL,EarthPath[i][1] );

    apcModel->SetPoint(2*i,sumLlow,EarthPathPREM[i][1] );

    apcModel->SetPoint(2*i+1,sumLup,EarthPathPREM[i][1]  );



  }

  sumL = 0;

  sumLup = 0; 

  std::cout << "Calculation" << EarthPathPREM.size() << " " << EarthPath.size()  <<std::endl;

  std::cout << "    " << std::endl;

  TGraph * TestModel1 = new TGraph(2*EarthPath.size());
  TestModel1->SetLineColor(kRed);
  TestModel1->SetMarkerStyle(45);
  TestModel1->SetMarkerColor(kRed);

  TGraph * TestModel2 = new TGraph();
  TestModel2->SetMarkerStyle(20);

  double sumGoal = 0;

  int kk = 0;

  for (int i = 0; i < EarthPath.size() ; ++i)
  {

    sumLup += EarthPath[i][0]; 
    sumLlow = sumLup - EarthPath[i][0];   


    //std::cout << " X             "<< sumLlow<< " "  << sumLmid << " " <<sumLup <<  std::endl;

    std::cout << " Interpolation "<< EarthPath[i][1] << " "  << EarthPathPREM[i][1]  <<  std::endl;


    TestModel1->SetPoint(2*i,sumLlow,100*(EarthPath[i][1]-EarthPathPREM[i][1])/EarthPathPREM[i][1] );

    TestModel1->SetPoint(2*i+1,sumLup,100*(EarthPath[i][1]-EarthPathPREM[i][1])/EarthPathPREM[i][1]);

    if ( 100*(EarthPath[i][1]-EarthPathPREM[i][1])/EarthPathPREM[i][1] > 0.01 )
    {
      sumGoal += (sumLup-sumLlow);

      TestModel2->SetPoint(2*kk,sumLlow,100*(EarthPath[i][1]-EarthPathPREM[i][1])/EarthPathPREM[i][1] );

      TestModel2->SetPoint(2*kk+1,sumLup,100*(EarthPath[i][1]-EarthPathPREM[i][1])/EarthPathPREM[i][1]); 

      ++kk;
    }

    //TestModel2->SetPoint(2*i,sumLlow,MyModel->Eval(sumLlow)  );

    //TestModel2->SetPoint(2*i+1,sumLup,MyModel->Eval(sumLup)  );

  }




  OscProbPath.close();
  ModelPath.close();


  ThePaths->Add(MyModel);
  ThePaths->Add(apcModel);
  ThePaths->Add(TestModel1);
  ThePaths->Add(TestModel2);

 std::cout << " Check here " << sumGoal << " " << test.Baseline_LLVP <<  std::endl;


 //Histogram showig distance trough LLVP







  /*



    for (int i = 0; i < paths.size() ; ++i)
  {

    sumLup += paths[i].length ;
    sumLlow = sumLup - paths[i].length ;   
    std::cout << i << " " << sumL << " " << paths[i].length  << " " <<  paths[i].density  << " " << paths[i].zoa << std::endl;

    OscProbPath << sumL << " " << paths[i].length << " " << paths[i].density << std::endl;

//    apcModel->SetPoint(i,sumL,paths[i].density  );

    apcModel->SetPoint(2*i,sumLlow,paths[i].density  );

    apcModel->SetPoint(2*i+1,sumLup,paths[i].density  );




  }







  std::string NuTomoPath= "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation";
  
  std::string DistanceFolder = "/SimulationResults/PreliminaryResults/BaselinesLLVPs/";

  std::string BaselineFile = "BaselineLLVP_"+llvpShape+"_"+std::to_string(zenithBins)+"zenithBins"+"_"+ std::to_string(static_cast<int>(zenithMin)) +"-"+ 
                                 std::to_string(static_cast<int>(zenithMax))+"Zen_"+ std::to_string(azimuthBins)+ std::to_string(static_cast<int>(azimuthMin)) +"-"+ 
                                 std::to_string(static_cast<int>(azimuthMax))+"Az.csv";



  TH2D* LLVPBaseLineHIst = new TH2D("LLVPBaseLine", "LLVP Base Line", zenithBins, cosZenithMin, cosZenithMax, azimuthBins, azimuthMin, azimuthMax);

    for (int i = 1; i <= LLVPBaseLineHIst->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= LLVPBaseLineHIst->GetNbinsY(); ++j)
        {
            double cthi = LLVPBaseLineHIst->GetXaxis()->GetBinCenter(i);
            double phij = LLVPBaseLineHIst->GetYaxis()->GetBinCenter(j)*TMath::Pi()/180.0;

            Baselines.SetModel(path2prem);

            Baselines.SetDirection(cthi, phij);
  
            Baselines.SetPile( true, llvpShape, 3.0, 0.0);


            std::vector<std::vector<double>> InLLVP_Path = Baselines.Create3DPath( );
          
            double  Lij= Baselines.Baseline_LLVP;

            LLVPBaseLineHIst->SetBinContent(i, j, Lij);

        }
    }

  ExportToCSV(LLVPBaseLineHIst, NuTomoPath+DistanceFolder+BaselineFile);


  




 
 


  TCanvas * c2 = new TCanvas();

  LLVPBaseLineHIst->Draw("COLZ");

  LLVPBaseLineHIst->SetStats(0);

  c2->Modified(); c2->Update();
  

  
  TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();
  rc2->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

  */

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


void ExportToCSV(TH2D* hist, std::string filename)
{


    std::ofstream outfile(filename);

    // Recorre los bins y guarda el contenido


    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        for (int j = 1; j <= hist->GetNbinsY(); ++j)
        {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);

            double N = hist->GetBinContent(i, j);
     

            outfile << x << "," << y << "," << N << std::endl;
        }
    }

    outfile.close();
}
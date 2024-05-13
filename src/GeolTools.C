/*
In this code:

PremMatrix = [R, RHO, Z/A, LAYER] 
PathMatrix = [L, RHO, Z/A, LAYER] 

R --- Column vector of layer radius
L --- Column vector of Neutrino Baselines


*/

#include "GeolTools.h"

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
#include "TMath.h"
#include <Math/SVector.h>
#include "TCanvas.h"
#include "TSystem.h"

//Geometry manager

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoVolume.h"
#include "TGLViewer.h"

using namespace std;



//#include <vector> 
std::vector< std::vector<double> > SortPREMdata( std::string PREM_MODEL = "prem_15layers.txt" )
{
   std::string PREM_PATH = "/home/dehy0499/OscProb/PremTables/"+PREM_MODEL;

   //--- open PREM model 
   std::vector< std::vector<double> > PremMatrix; // Matrix data from Prem model
   std::vector<double> PremData(4);               // [R, rho, z/a, layer]

   // Variables for storing table rows
   double radius, density, zoa, layer;

   std::ifstream PREM_DATA;

   PREM_DATA.open(PREM_PATH);

      // Loop over table rows
   while(PREM_DATA >> radius >> density >> zoa >> layer)
   {

      PremData[0] = radius;
      PremData[1] = density;
      PremData[2] = zoa;
      PremData[3] = layer;

      PremMatrix.push_back(PremData);
    }

   return PremMatrix;

}


void PREM3D( std::vector< std::vector<double> > PremMatrix, bool LLVP, double h_llvp = 1000, double drho_llvp = 1.03, double LLVP_angle = 45.0 )
{

   
   int last = PremMatrix.size()-1; // Index of last Elment

   double Ro = PremMatrix.back()[0]; //Outermost Radius in the Models

   //--Define some transformation
   TGeoTranslation *tr1 = new TGeoTranslation("tr1", 0, 0, 1.0);

   TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);//Rotation
   
   TGeoCombiTrans *combi1 = new TGeoCombiTrans("combi1",0, 0, 1.0, rot1 );
   
   TCanvas *c1 = new TCanvas();
   
   TGLViewer *view = (TGLViewer*)gPad->GetViewer3D();
   
   TGeoManager *geom = new TGeoManager(); // Instance of object 

   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0); // World Material
 
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",9999, matVacuum); // World Medium
    
   TGeoVolume *world = geom -> MakeBox("TOP", Vacuum , 2, 2, 3);

   geom->SetTopVolume(world); //Make top medium the world
     

   //Inner part of Inner Core

   double radmin = 0.0;

   double radmax = PremMatrix[0][0]/Ro;

   double density = PremMatrix[0][1]; 
   
   double zoa = PremMatrix[0][2];
   
   double layer = PremMatrix[0][3];
  
  
    TGeoMaterial *matlay = new TGeoMaterial("matlay",1.0 ,zoa,density);
    //--- define some media
    TGeoMedium *medlay= new TGeoMedium("Layer Material",1.0, matlay);
 


   TGeoVolume *modlayer= geom -> MakeSphere("layer", medlay ,radmin,radmax);

   //modlayer->SetVisibility(kTrue);
   
   world->AddNode(modlayer,99999,tr1);
   
   for (int i = 1; i < PremMatrix.size(); ++i)
   {
      

      radmin = PremMatrix[i-1][0]/Ro;
      
      radmax = PremMatrix[i][0]/Ro;

      double radmid = (radmin + radmax)/2.0; 

      density = PremMatrix[i][1]; 
      
      zoa = PremMatrix[i][2];
      
      layer = PremMatrix[i][3];

      TGeoMaterial *matlay = new TGeoMaterial("mat",1.0 ,zoa,density);
       
      TGeoMedium *medlay = new TGeoMedium("layer Material",layer, matlay);

      TGeoVolume *modlayer = geom -> MakeSphere("layer", medlay , radmin, radmax);

      //std::cout << " radmin " << PremMatrix[i-1][0] << " radmax " << PremMatrix[i][0] << " " << radmid << " R " << Ro <<std::endl;

      modlayer->SetLineColor(kOrange);

     // modlayer->SetVisibility(kTrue);
      
      // Insert LLVPs Sectors

      //int j = 1; //LLVP_index

       if ( LLVP )
       { 

          std::cout << " LLVP IS ACTIVATED " << std::endl;
          
          double local_rho = drho_llvp*density; //Fixed density

          double llvp_rmax = 3480.0 + h_llvp; //Fixed density

          if (radmin*Ro >= 3480.0 && radmax*Ro <= llvp_rmax)
          {

          std::cout << "****************************-*-*-*-*-*" << radmin*Ro << " " << radmax*Ro << " " << (radmax-radmin)*Ro << " " << h_llvp/3 <<  std::endl;
          
             
             TGeoMaterial *matblob= new TGeoMaterial("matblob",1.0 ,zoa,local_rho);
                
             TGeoMedium *medblob = new TGeoMedium("LLVPs Material",7, matblob);
               
             TGeoVolume *LLVP = geom -> MakeSphere("LLVP", medblob , radmin, radmax , 0,LLVP_angle,0,360 );
                
             modlayer->AddNode(LLVP,i+1,rot1);

          }

       }


       world->AddNode(modlayer,i+1,tr1);
   
   }

   geom->CloseGeometry();

   //   world->SetLineColor(kBlack);

   //   gGeoManager->SetTopVisible(); // the TOP is invisible

   //   world->Draw();

}


std::vector< std::vector<double> >  NuPATHS3D (std::string PREM_MODEL, double zen, double azi,  bool LLVP )
{
   gSystem->Load("libGeom");

   std::vector< std::vector<double> > PathMatrix;
  
   std::vector<double> PathData(4); // Row [L, rho, z/a, layer]

  
   std::vector< std::vector<double> > PremMatrix = SortPREMdata(PREM_MODEL); // Sort PREM model into a readable matrix

   // Displaying the PREM Model
   
   /*
    for (int i = 0; i < PremMatrix.size(); i++) 
    { 

      cout << PremMatrix[i][0] << " " << PremMatrix[i][1] << " "<< PremMatrix[i][2] << " "<< PremMatrix[i][3] << " "; 
      cout << endl; 
    
    } 
    
    std::cout << std::endl; 

    */



   PREM3D( PremMatrix, LLVP );

   int last = PremMatrix.size()-1; // Index of last Elment

   double Ro = PremMatrix.back()[0]; //Outermost Radius in the Models


   //Calculate Paths inside the Earth

   //Direction of neutrino in spherical coordiates: https://mathworld.wolfram.com/SphericalCoordinates.html
   //double zen = 180.0; // zen (90,180]

   double pol = zen; // 0 to Pi
   
   //double azi = 0.0; // 0 to 2*Pi
   
   //Initial Position
   double x = sin(pol)*cos(azi);
   double y = sin(pol)*sin(azi);
   double z = cos(pol);

   double n = sqrt(x*x + y*y + z*z);

   //Unit vector
   double nx = x/n;
   double ny = y/n;
   double nz = z/n;

   //Origin of the Sphere

   double zo = Ro/Ro;

   //Initial Point must be in the outermost layer -> Intersection Line spehere

   ROOT::Math::SVector<double, 3> u(nx,ny,nz); // Direction of neutrino
   ROOT::Math::SVector<double, 3> o(0.0,0.0,0.0); // Origin of the line
   ROOT::Math::SVector<double, 3> cc(0,0, zo); // Origin of the Sphere 





   double r = Ro/Ro;

   double a = ROOT::Math::Dot( u, u);
   double b = 2*ROOT::Math::Dot( u, o - cc);
   double c = ROOT::Math::Dot( o-cc, o-cc) - r*r ;

   double d1 = (-b+sqrt(b*b -4*a*c))/(2*a);
   double d2 = (-b-sqrt(b*b -4*a*c))/(2*a);

   //std::cout << " d1 " << d1 << " d2 " << d2 << std::endl; 

   //std::cout << "Initial direction " << -nx << " " << -ny << " " << -nz << std::endl;


    gGeoManager->InitTrack(d1*nx,d1*ny,d1*nz, -u[0] , -u[1], -u[2]); // d*u is the point in the sphere

    // TPolyLine3D *l = new TPolyLine3D();

    int i = 0;

    double sumL = 0;

    double Ltot = abs(2*Ro*cos(zen*TMath::Pi()/180.0));

    while (!gGeoManager->IsOutside ())
   {
      
         

         int nodeID= gGeoManager->GetCurrentNode()->GetIndex();

         const Double_t *cpoint = gGeoManager->GetCurrentPoint();

         const char *path = gGeoManager->GetPath();

         TGeoVolume *cvol = gGeoManager->GetCurrentVolume();

         TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial(); //Material of the current Boundary
            
         gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.

         Double_t Li = gGeoManager->GetStep(); //Baseline Segement

         Double_t safety = gGeoManager->GetSafeDistance(); //Distance to the next boundary/Neutrino Baseline

         if (cvol->GetMedium()->GetId() >= 10)
         {
            continue;
         }

         else {

         PathData[0] = Li*Ro; //To the next boundary

         PathData[1] = cmat->GetDensity();
         
         PathData[2] = cmat->GetZ();
         
         PathData[3] = cvol->GetMedium()->GetId(); // To Be adjusted!!!

          }

         sumL += Li*Ro;

         
      if (sumL <= Ltot) 
      {
         std::cout << "Inside Of Geometry" << std::endl;

         std::cout << "Current path is: " << path << " r "  << sqrt(cpoint[0]*cpoint[0]+ cpoint[1]*cpoint[1] + cpoint[2]*cpoint[2])<< std::endl;

         std::cout << "L " << PathData[0]  << " rho "  << PathData[1] << " zoa " << PathData[2] << " sumL " << sumL << " | " << Ltot <<std::endl;

         std::cout << std::endl; 
         

         // l->SetPoint( i , cpoint[0],cpoint[1],cpoint[2]);

         i += 1;

      }

      else
      {

         std::cout << "Outside Of Geometry" << std::endl;

      }
      
         PathMatrix.push_back(PathData);
      
   } 

   std::cout << std::endl;

   // Displaying Path data

   std::cout<< "Zenith " << zen << " Azimuthal " << azi*180.0/TMath::Pi() << std::endl;

   std::cout << std::endl;
   
    for (int i = 0; i < PathMatrix.size(); i++) { 
        //for (int j = 0; j < PremMatrix[i].size(); j++) 
            std::cout << PathMatrix[i][0] << " "<< PathMatrix[i][1] << " "<< PathMatrix[i][2] << " "<< PathMatrix[i][3] << " "; 
        std::cout << std::endl; 
    } 
    
    std::cout << std::endl; 
 /*
   l->SetLineColor(6);
   l->SetLineWidth(3);
   l->Draw("same");
  */ 
   return PathMatrix;
 }



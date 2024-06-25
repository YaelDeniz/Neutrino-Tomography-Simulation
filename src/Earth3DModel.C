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




//Cern Root
#include <math.h>
#include "TMath.h"
#include <Math/SVector.h>
#include "TCanvas.h"
#include "TSystem.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TView.h"


//Geometry manager

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoVolume.h"
#include "TGLViewer.h"


using namespace std;

//Set the Prem Model

void Earth3DModel::SetModel(std::string model)
{
    filename = model;
} 
  
//Set the Detector Location

void Earth3DModel::SetDetector(double Position[3]) 
{
   Det[0] = Position[0];
   Det[1] = Position[1];
   Det[2] = Position[2];
}

//Set Neutrino Direction

void Earth3DModel::SetDirection(double cth , double phi ) //Input in radias -> trasnformed back to degrees
{
    //zenith =  acos(cth); /*Polar direction*/

      // Check if the total number of elements matches
    if (cth < -1 || cth > 1) {
        throw std::invalid_argument("Invalid argument for neutrino direction, zenith direction must be in terms of cos(zen)");
    }

    else
    {

    zenith = TMath::Pi() - acos(cth); /*Polar direction*/


    // double th = TMath::Pi()*( 1-(zen/180.0) ); //Polar angle in spherical coordinates in rads [pi/2, pi] 

    
    azimuth = phi*TMath::Pi()/180.0;   /*Azimuthal direction in degrees*/
    
    std::cout <<"Direction: " <<  zenith*180.0/TMath::Pi() << " " << azimuth*180.0/TMath::Pi() << std::endl;

    }
} 

//Define if LLVP exist

void Earth3DModel::ActiveHeterogeneity( bool value , std::string shape) 
{ 
  Anomaly = value; 
  AnomalyShape = shape;
} // Activate the LLVPs


// Create a matrix from PremTables
std::vector< std::vector<double> > Earth3DModel::GetPremData( std::string PREM_MODEL  )
{
   std::string PREM_PATH = "/home/dehy0499/OscProb/PremTables/"+PREM_MODEL;

   //--- open PREM model 
   std::vector< std::vector<double> > PremMatrix; // Matrix data  form of Prem tables

   std::vector<double> PremRow(4);               // [radius density z/a, layer ID]

   double radius, density, zoa, layer; // Variables for storing Prem table values

   std::ifstream PREM_DATA;

   PREM_DATA.open(PREM_PATH);

      // Loop over table rows
   while(PREM_DATA >> radius >> density >> zoa >> layer)
   {

      PremRow[0] = radius;
      PremRow[1] = density;
      PremRow[2] = zoa;
      PremRow[3] = layer;

      PremMatrix.push_back(PremRow);
    }

   return PremMatrix;
}

//Label PREM Layers

int Earth3DModel::LabelLayer (double radius)
{
  int id;

  if( 0 <= radius && radius <= 1221.5) { id = 0;}
  else if(1221.5 < radius && radius <= 3480.0) {id = 1;}
  else if(3480.0 < radius && radius <= 5701.0) {id = 2;}
  else if(5701.0 < radius && radius <= 6346.6) {id = 3;}
  else if(6346.6 < radius && radius <= 6368.0) {id = 4;}
  else if(6368.0 < radius && radius <= 6371.0) {id = 5;}
  else if(6371.0 < radius && radius <= 6390.0) {id = 6;}
  // else {id = 0;}
  return id;
}


//CONSTRUCT LLVP

void Earth3DModel::CreatePanCake(std::vector<TGeoVolume*> LAYER, std::vector< std::vector<double> > PremMatrix, std::vector< std::vector<double> > check )
{

  //Baker----------------------------------------------------------------------------------------------------------------------
  
  std::vector<TGeoMaterial*> LLVPMat; //Vector to define material for each LLVP segment ( A X% anomaly in local Density or Z/A)

  std::vector<TGeoMedium*> LLVPMed;   //Vector to define medium for each LLVP segment

  std::vector<TGeoVolume*> LLVPLAYER; //Vector to define Volume for each LLVP segment
  
  TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);// Some Geometrical trasformation that move LLVPs to the right place

  double PileThickness = 1000;

  double Bottom,Top;

  double Localdensity;

  double Localchem;

  Bottom = 3480;
  
  Top = 3480 + PileThickness;

   // Segment 1 test
   double slabLength =0;

   double InnerR = 0;
   double OuterR = 0;

   int PileId = 0;

   std::cout << " " << std::endl;
   std::cout << "Pile radius: " << Bottom << " - " << Top << std::endl;

   double lowerSegment = Bottom; //Lower Radius of a segment in Pile section

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInPile = PremMatrix[i][0]; //Lower boundary of current layer

      if (RadiusInPile > Bottom )
      {

          char integer_string[32];
            
          sprintf(integer_string, "%d", i);//Layer Index

          char llvpmaterial_string[64]="LLVPMaterial";
          strcat(llvpmaterial_string, integer_string);
          const char *llvpMatName = llvpmaterial_string;
          

          char llvpmedium_string[64]="LLVPMedium";
          strcat(llvpmedium_string, integer_string);
          const char *llvpMedName = llvpmedium_string;
            
          char llvplayer_string[64]="LLVPLayer";
          strcat(llvplayer_string, integer_string);
          const char *llvpLayerName= llvplayer_string;
        
         slabLength = slabLength + (RadiusInPile - lowerSegment ); //Calculate segement lenght 

         if (slabLength > PileThickness)
         { 

            double excess = PileThickness - slabLength;
            InnerR= lowerSegment;
            OuterR= RadiusInPile + excess;
  
          //std::cout << "top of lower segment: " << InnerR << " " << OuterR << std::endl;
 

            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
            Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(llvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(llvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(llvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            LLVPLAYER[PileId]->SetLineColor(kBlue);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            break;
         }
    
         else
         {
            InnerR = lowerSegment;
            OuterR = RadiusInPile;
            
            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
          
            //std::cout << "Radius of lower section's segements: "<< InnerR << " " << OuterR << std::endl;

            Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(llvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(llvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(llvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            LLVPLAYER[PileId]->SetLineColor(kBlue);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            ++PileId;

         }

         lowerSegment = RadiusInPile;
      }
   }
}

void Earth3DModel::CreateCake(std::vector<TGeoVolume*> LAYER, std::vector< std::vector<double> > PremMatrix, std::vector<std::vector<double>> check )
{

  //Baker----------------------------------------------------------------------------------------------------------------------
  
  std::vector<TGeoMaterial*> LLVPMat; //Vector to define material for each LLVP segment ( A X% anomaly in local Density or Z/A)

  std::vector<TGeoMedium*> LLVPMed;   //Vector to define medium for each LLVP segment

  std::vector<TGeoVolume*> LLVPLAYER; //Vector to define Volume for each LLVP segment
  
  //Naming Convention

  //Generate segment names
  //char integer_string[32];
    
  //sprintf(integer_string, "%d", i+1);

  //char llvpmaterial_string[64]="LLVPMaterial";
  //char llvpmaterial_string[64];
  //strcat(llvpmaterial_string, integer_string);
  //const char *llvpMatName = llvpmaterial_string;
  //const char *llvpMatName;

  //char llvpmedium_string[64]="LLVPMedium";
  //strcat(llvpmedium_string, integer_string);
  //const char *llvpMedName = llvpmedium_string;
  //const char *llvpMedName;
    
  //llvplayer_string[64]="LLVPLayer";
  //strcat(llvplayer_string, integer_string);
  //const char *llvpLayerName= llvplayer_string;
  //const char *llvpLayerName;


   TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);// Some Geometrical trasformation that move LLVPs to the right place

   double PileThickness = 1000;

   double thicknessOfSection = PileThickness/3.0; // Salabas have same thickness

   double lowerBottom, lowerTop, midBottom, midTop, upperBottom, upperTop;

   double Localdensity;

   double Localchem;

   lowerBottom = 3480;

   lowerTop = 3480 + thicknessOfSection;

   midTop = 3480 + 2*thicknessOfSection;

   upperTop = 3480 + 3*thicknessOfSection;

   // Segment 1 test
   double slabLength =0;

   double InnerR = 0;
   double OuterR = 0;

   int PileId = 0;

   std::cout << " " << std::endl;
   std::cout << " Segement 1 radius: " << lowerBottom << " - " << lowerTop << std::endl;

   double lowerSegment = lowerBottom; //Lower Radius of a segment in Pile section

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInBottom = PremMatrix[i][0]; //Lower boundary of current layer

      if (RadiusInBottom > lowerBottom )
      {

          char lowinteger_string[32];
            
          sprintf(lowinteger_string, "%d", i);//Layer Index

          char lowllvpmaterial_string[64]="lowLLVPMaterial";
          strcat(lowllvpmaterial_string, lowinteger_string);
          const char *lowllvpMatName = lowllvpmaterial_string;
          

          char lowllvpmedium_string[64]="lowLLVPMedium";
          strcat(lowllvpmedium_string, lowinteger_string);
          const char *lowllvpMedName = lowllvpmedium_string;
            
          char lowllvplayer_string[64]="lowLLVPLayer";
          strcat(lowllvplayer_string, lowinteger_string);
          const char *lowllvpLayerName= lowllvplayer_string;
        
         slabLength = slabLength + (RadiusInBottom - lowerSegment ); //Calculate segement lenght 

         if (slabLength > thicknessOfSection)
         { 

            double excess = thicknessOfSection - slabLength;
            InnerR= lowerSegment;
            OuterR= RadiusInBottom + excess;
  
          //std::cout << "top of lower segment: " << InnerR << " " << OuterR << std::endl;
 
            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
            
            Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(lowllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(lowllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(lowllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            LLVPLAYER[PileId]->SetLineColor(kBlue);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            midBottom = OuterR;
            
            break;
         }
    
         else
         {
            InnerR = lowerSegment;
            OuterR = RadiusInBottom;
            
            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
          
            //std::cout << "Radius of lower section's segements: "<< InnerR << " " << OuterR << std::endl;

            Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(lowllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(lowllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(lowllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, aWidth ,0 ,360 ) ); // DEFINE LLVP Segment

            LLVPLAYER[PileId]->SetLineColor(kBlue);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            ++PileId;

         }

         lowerSegment = RadiusInBottom;
      }
   }

  
   //Segment 2------------------------------------------------------------------------------------------------------------------------
   PileId=0;

   slabLength =0;
   InnerR = 0;
   OuterR = 0;
   
   std::cout << " " << std::endl;
   std::cout << " Segement 2 radius: " << midBottom << " - " << midTop << std::endl;
   
   double midSegment = midBottom; //Lower Radius of a segment in Pile section

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInMiddle = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (RadiusInMiddle > midBottom )
      {
        char midinteger_string[32];
            
        sprintf(midinteger_string, "%d", i);//Layer Index

        char midllvpmaterial_string[64]="midLLVPMaterial";
        strcat(midllvpmaterial_string, midinteger_string);
        const char *midllvpMatName = midllvpmaterial_string;
          

        char midllvpmedium_string[64]="midLLVPMedium";
        strcat(midllvpmedium_string, midinteger_string);
        const char *midllvpMedName = midllvpmedium_string;
            
        char midllvplayer_string[64]="midLLVPLayer";
        strcat(midllvplayer_string, midinteger_string);
        const char *midllvpLayerName= midllvplayer_string;

         std::cout << " Radius of layer above middle bottom  " << RadiusInMiddle << std::endl;


         slabLength = slabLength + (RadiusInMiddle - midSegment ); //Calculate segment lenght 

         if (slabLength > thicknessOfSection)
         {
            double midExcess = thicknessOfSection-slabLength;
            InnerR= midSegment;
            OuterR= RadiusInMiddle + midExcess;


            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
            
            
            Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(midllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(midllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            double midwidth= (2.0/3.0)*aWidth;

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(midllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, midwidth ,0 ,360 ) ); // DEFINE LLVP Segment

            LLVPLAYER[PileId]->SetLineColor(kGreen);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            upperBottom = OuterR;

            break;
         }

         else
         {
          InnerR = midSegment;
          OuterR = RadiusInMiddle;

          double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
          double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
          std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
          std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
          
          Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

          Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
  
          LLVPMat.push_back( new TGeoMaterial(midllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material

          LLVPMed.push_back( new TGeoMedium(midllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

          double midwidth= (2.0/3.0)*aWidth;

          LLVPLAYER.push_back( gGeoManager -> MakeSphere(midllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, midwidth ,0 ,360 ) ); // DEFINE LLVP Segment

          LLVPLAYER[PileId]->SetLineColor(kGreen);

          LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1); 

          ++PileId; 

         }
         midSegment = RadiusInMiddle;
      }

   }
   //Segment 3------------------------------------------------------------------------------------------------------------------------
   PileId=0;
   slabLength =0;
   InnerR = 0;
   OuterR = 0;
   
   std::cout << " " << std::endl;
   
   std::cout << " Segement 3 radius: " << upperBottom << " - " << upperTop << std::endl;

   double upperSegment = upperBottom;

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInTop = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (RadiusInTop > upperBottom )
      {


        char upinteger_string[32];
            
        sprintf(upinteger_string, "%d", i);//Layer Index

        char upllvpmaterial_string[64]="upLLVPMaterial";
        strcat(upllvpmaterial_string, upinteger_string);
        const char *upllvpMatName = upllvpmaterial_string;
          

        char upllvpmedium_string[64]="upLLVPMedium";
        strcat(upllvpmedium_string, upinteger_string);
        const char *upllvpMedName = upllvpmedium_string;
            
        char upllvplayer_string[64]="upLLVPLayer";
        strcat(upllvplayer_string, upinteger_string);
        const char *upllvpLayerName= upllvplayer_string;

         std::cout << " Radius of layer above uppper bottom  " << RadiusInTop << std::endl;

         slabLength = slabLength + (RadiusInTop - upperSegment );

         //r.push_back(dh-mslabL);

         if (slabLength > thicknessOfSection)
         {
            double topExcess = thicknessOfSection-slabLength;
            InnerR= upperSegment;
            OuterR= RadiusInTop + topExcess;

            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
            
            
            Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(upllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(upllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            double upwidth= (1.0/3.0)*aWidth;

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(upllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, upwidth ,0 ,360 ) ); // DEFINE LLVP Segment

            LLVPLAYER[PileId]->SetLineColor(kRed);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            LLVPMat.clear();
            
            LLVPMed.clear();
            
            LLVPLAYER.clear();

            break;
         }

         else
         {
            InnerR = upperSegment;
            OuterR = RadiusInTop;

            double rhotest = LAYER[i]->GetMedium()->GetMaterial()->GetDensity();
            double zoatest = LAYER[i]->GetMedium()->GetMaterial()->GetZ();
            std::cout <<"EARTH| rmin: " << check[i][0] << " rmax: " << check[i][1] <<" rho:" << rhotest << " zoa: " << zoatest << std::endl;
            std::cout<<"LLVP| rmin: " <<InnerR << " rmax: " << OuterR<< " rho:" << PremMatrix[i][1] << " zoa: " << PremMatrix[i][2]<<std::endl;
            
            
            
            Localdensity = (1.0 + drho/100.0)*PremMatrix[i][1];

            Localchem = (1.0 + dzoa/100.0 )*PremMatrix[i][2];
    
            LLVPMat.push_back( new TGeoMaterial(upllvpMatName, 1, Localchem, Localdensity) ); // Create LLVP material
  
            LLVPMed.push_back( new TGeoMedium(upllvpMedName,1,LLVPMat[PileId]) ); // Create LLVP medium

            double upwidth= (1.0/3.0)*aWidth;

            LLVPLAYER.push_back( gGeoManager -> MakeSphere(upllvpLayerName, LLVPMed[PileId] , InnerR , OuterR , 0, upwidth ,0 ,360 ) ); // DEFINE LLVP Segment

            LLVPLAYER[PileId]->SetLineColor(kRed);

            LAYER[i]->AddNode(LLVPLAYER[PileId],1,rot1);  

            ++PileId;

         }
         upperSegment = RadiusInTop;
      }

   }
   
}






// CONSTRUCT THE 3D MODEL

std::vector<std::vector<double>> Earth3DModel::Earth3DPath( double th, double phi, std::string MODEL )
{
  //Read Prem Data

  gSystem->Load("libGeom");

  std::cout << "Earth Model: " << MODEL << std::endl;

  std::vector< std::vector<double> > PremMatrix = GetPremData(MODEL); // Sort PREM model into a readable matrix

  std::vector< std::vector<double> > LLVPMatrix = GetPremData(MODEL); // A Copy of PREM model matrix to construct LLVPs

  double rPREM = PremMatrix.back()[0]; // Outermost Layer in the PREM model

  std::vector< std::vector<double> > EarthPath;  // Store Paths inside the Earth for oscillation calculations

  std::vector< int > EarthPathId;  // Store ID layer for each Paths inside the Earth for oscillation calculations

  EarthCanvas = new TCanvas("3D Earth", "3D Earth",0,0,600,600);
 
  new TGeoManager("EarthModel3D", "Simple 3D geometry");

  //Define world (Volume container)

  //Defining Vacuum Medium
  TGeoMaterial *VAC = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
  
  TGeoMedium *vac = new TGeoMedium("VACUUM",99,VAC);
  
  TGeoVolume *top = gGeoManager->MakeBox("TOP",vac,rPREM+5000,rPREM+5000,rPREM+5000);

  gGeoManager->SetTopVolume(top);

  //Create Earth

  /* Each Layer is represented of a Spherical Shell defined by a Inner Radius (rmin) and Outer Radius(rmax).
  Each Layer Medium is defined by a Material, a Material is created from a specific density value and Z[proton]/A[neutron + proton] (zoa)*/
  double rmin, rmax, density, zoa, layer; 

  std::vector<TGeoMaterial*> MAT; //Vector to define material for each layer en in PREM (Density & Z/A)
  std::vector<TGeoMedium*>  MED;  //Vector to define medium for each layer en in PREM 
  std::vector<TGeoVolume*> LAYER; //Vector to define a  volume for each layer en in PREM
  std::vector<std::vector<double>> check;

    std::cout << "***************************THE EARTH***************************" << std::endl; 


  for (int i = 0; i < PremMatrix.size(); ++i)
   {
      if (i==0) { rmin = 0; }

      else { rmin = PremMatrix[i-1][0]; }

      rmax = PremMatrix[i][0]; //Layer Outer radius from PREM 

      density = PremMatrix[i][1];  //Layer density from PREM
      
      zoa = PremMatrix[i][2]; //Layer Z/A from PREM
      
      layer = PremMatrix[i][3]; //Layer ID number from PREM

      //Generate Names for Material, Medium and Layer
      char integer_string[32];

      sprintf(integer_string, "%d", i+1);

      char material_string[64]="Material";
      strcat(material_string, integer_string);
      const char *MatName = material_string;

      char medium_string[64]="Medium";
      strcat(medium_string, integer_string);
      const char *MedName = medium_string;

      char layer_string[64]="Layer";
      strcat(layer_string, integer_string);
      const char *LayerName= layer_string;

      //Definiton of world "Daugther" volumes (i.e., Earth Layers)

      MAT.push_back( new TGeoMaterial( MatName, 1.0 , zoa , density) ); //Define Material for ith Layer 
       
      MED.push_back( new TGeoMedium( MedName , 1 , MAT[i]) ); //Define Medium for ith Layer 

      check.push_back({rmin,rmax});

      LAYER.push_back( gGeoManager->MakeSphere(LayerName,MED[i], rmin,rmax,0,180,0,360) ); // Define Volume for the ith Layer

      if (rmax == 3480.0 ) { LAYER[i]->SetVisibility(kTRUE);} //Outer Core is visible

      else { LAYER[i]->SetVisibility(kFALSE); }
  }

      LAYER[LAYER.size()-1]-> SetVisibility(kTRUE); // Crust is vissible

  /*

  // SECTION TO BE DELETED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double R_cmb = 3480.0;
  double zentest = 180-TMath::ASin( (R_cmb)/6368.0 )*(180.0/TMath::Pi()) ; // max 180
  if (zen  > zentest )
  {
    Anomaly = false;
  }
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  */

  std::cout << "***************************THE LLVP***************************" << std::endl; 
//-+++++++++++++++++++++++++++++++++++++++++++++-----++-+-+-+-+--+

  if (Anomaly)
  {

  if (AnomalyShape == "cake")
  {
     CreateCake(LAYER, LLVPMatrix,check );
  }

  else if (AnomalyShape == "pancake")
  {
     CreatePanCake(LAYER, LLVPMatrix, check );
  }
  

  }

  //-+++++++++++++++++++++++++++++++++++++++++++++-----++-+-+-+-+--+
  //Add Earth Layers to TOP Volume

  for (int i = 0; i < PremMatrix.size(); ++i)
  {
        top->AddNode(LAYER[i],1);
  }

  gGeoManager->CloseGeometry(); // Finish Geometry

  
  
  gGeoManager->SetTopVisible(); // the TOP is invisible

  
  
 
  top->Draw();

  TView *view = gPad->GetView();
  
  view->ShowAxis();

  
 
  // Tracking----------------------------------------------------------------------------------------------------------------------------------------

  //GetInitialPosition(zen,azi)
  //Calculate Paths inside the Earth

   //Direction of neutrino in spherical coordiates: https://mathworld.wolfram.com/SphericalCoordinates.html
 // double zen = 180.0; // zenith angle respect to the detector location in degrees (90,180]

 // double azi = 0.0;  //Azimuthal angle respect to the detector location in degrees [0,360]

  
  // double th = TMath::Pi()*( 1-(zen/180.0) ); //Polar angle in spherical coordinates in rads [pi/2, pi] 

  //double phi = TMath::Pi()*(azi)/180.0; //Azimuthal angle in spehericla coordinates in rads [0, 2*pi]
   


  double r = rPREM; // Max Radius

  //double Det[3]= {0.0,0.0,-6371.0}; //Detector location in cartesian coordinates
  
  //Incoming neutrino direction Cosines (respect to detector location)
  double nx = sin(th)*cos(phi);
  double ny = sin(th)*sin(phi);
  double nz = cos(th);


  ROOT::Math::SVector<double, 3> u(nx, ny, nz); // Direction of incoming neutrino
  ROOT::Math::SVector<double, 3> o(Det[0], Det[1], Det[2]); // Origin of the line (Set to be Detector location)
  ROOT::Math::SVector<double, 3> cc(0,0,0); // Origin of the Sphere 
 
  /* To calculate the initial position of the neutrino we need to compute 
  the intersection betweem 3D line (Neutrino direction) and Sphere(Earth) */

  double a = ROOT::Math::Dot( u, u);
  double b = 2*ROOT::Math::Dot( u, o - cc);
  double c = ROOT::Math::Dot( o-cc, o-cc) - r*r ;

  //Initial position of the neutrino... It should be located somewhere in the Atmosphere (Outer Layer of the Earth)
  double d1 = (-b+sqrt(b*b -4*a*c))/(2*a); // Distance 1 from the origin of the Line (Detector Location)
  double d2 = (-b-sqrt(b*b -4*a*c))/(2*a); // Distance 2 from the origin of the Line (Detector Location)

  
  //Coordinates of the Neutrino Positionin the Earth

  double xo = o[0]+d1*u[0];
  double yo = o[1]+d1*u[1];
  double zo = o[2]+d1*u[2];

  // Neutrino Direction

  double ni = -u[0];
  double nj = -u[1];
  double nk = -u[2];
  
  gGeoManager->InitTrack(xo, yo, zo, ni, nj, nk); // d*u is the point in the sphere

  TPolyMarker3D *l1 = new TPolyMarker3D(2,2); //Markers indicating DetLoc and Inical Neutrino
  l1->SetPoint( 0 , o[0],o[1],o[2]); //Detector
  l1->SetPoint( 1 ,xo,yo,zo); //Incoming Neutrino



  std::cout << o[0] << " " << o[1] << " " << o[2] << std::endl;
  std::cout << xo << " " << yo << " " << zo << "R: " << sqrt(xo*xo + yo*yo + zo*zo) <<std::endl;

  TPolyLine3D *l2 = new TPolyLine3D(); // Lines that represent neutrino Paths.
  


  int i = 0;

  double Dx2, Dy2, Dz2, Dnorm;

  double sumL = 0; //Track Total Baseile

  //Neutrino Propagation inside the Earth
  
  std::cout << "***************************Setting Neutrino***************************" << std::endl; 
  double v = 0;
  while (!gGeoManager->IsOutside ())
  {  
       ++v;
       std::cout << v << std::endl;

        // std::cout << "Starting loop" << std::endl;
         
         //gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
      
         int nodeID= gGeoManager->GetCurrentNode()->GetIndex();

         const Double_t *cpoint = gGeoManager->GetCurrentPoint(); // Current nuetrino Position

         Dx2 = (cpoint[0]-Det[0])*(cpoint[0]-Det[0]);
         Dy2 = (cpoint[1]-Det[1])*(cpoint[1]-Det[1]);
         Dy2 = (cpoint[2]-Det[2])*(cpoint[2]-Det[2]);

         Dnorm = round(sqrt( Dx2 + Dy2 + Dz2 ));


         if (Dnorm <= 0.0 )
         {

          //std::cout << " At detector" << std::endl;
          break;
         
         }
         

         const char *path = gGeoManager->GetPath();

         TGeoVolume *cvol = gGeoManager->GetCurrentVolume();

         TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial(); //Material of the current Boundary


            

         //Double_t Li = gGeoManager->GetStep(); //Baseline Segement

         Double_t safety = gGeoManager->GetSafeDistance(); //Distance to the next boundary/Neutrino Baseline

         gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
          
         Double_t Li = gGeoManager->GetStep(); //Baseline Segement

         double R_i = sqrt(cpoint[0]*cpoint[0]+ cpoint[1]*cpoint[1] + cpoint[2]*cpoint[2]); //Current Radius

        if (cvol->GetMedium()->GetId() == 99)
        {

          //std::cout << " OUTSIDE ************************************ " << std::endl;
          continue;

        }


         EarthPath.push_back({Li, cmat->GetDensity(), cmat->GetZ(),R_i}); //Store Path information [Baseline segment, Density, Z/A, R_i]
        
         EarthPathId.push_back(LabelLayer(rPREM - sumL)); // Store layer ID

         sumL = sumL + Li; //Track Total Segment (Max = -2*R_earth*cos(zen))

        l2->SetPoint( i , cpoint[0],cpoint[1],cpoint[2]);




        //Path Infromation 

        

         //std::cout << "Current path is: " << path << " Layer radius: "  << R_i << " id: " << LabelLayer(R_i) <<  std::endl;

         //std::cout << "L_i " << Li  << " rho "  << cmat->GetDensity() << " zoa " << cmat->GetZ() << " medId " << cvol->GetMedium()->GetId() << std::endl;

         //std::cout << std::endl; 

         i += 1;

        
        //stop at detetor
        /*

         Dx2 = (cpoint[0]-o[0])*(cpoint[0]-o[0]);
         Dy2 = (cpoint[1]-o[1])*(cpoint[1]-o[1]);
         Dz2 = (cpoint[2]-o[2])*(cpoint[2]-o[2]);
         
         Dnorm = round(sqrt(Dx2 + Dy2 + Dz2));

         std::cout << Dnorm << " ******|***** End loop" <<std::endl;




        
         if ( Dnorm < 0.1  )
         {
           break;
         }
         */

         
  }

   l1->SetMarkerColor(6);
   l1->SetMarkerSize(3);
   l2->SetLineColor(14);
   l2->SetLineWidth(3);
   
   l1->Draw("same");
   l2->Draw("same");

   EarthCanvas->Modified();
   EarthCanvas->Update();


   if(Anomaly)
   {

    //EarthCanvas->Print("SimulationResults/3DEarthLLVP.png");

   } 


  //DrawPrem

  //TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);

     //Geometrical display
   //std::cout<< " " << std::endl;
   //std::cout << "Detector Geometrical settings  d1 " << d1 << " d2 " << d2 << std::endl;
   //std::cout << " Detector: "  << Det[0] << " " << Det[1] << " " << Det[2] << std::endl;
   //std::cout << " Neutrino: "  << xo << " " << yo << " " << zo << std::endl;
   //std::cout << " Neutrino: "  << sqrt(xo*xo +yo*yo + zo*zo) << std::endl;

   std::cout << "***************************world created***************************" << std::endl; 
   return EarthPath;


}

// Create LLVP

 std::vector<std::vector<double>> Earth3DModel::Create3DPath () { return Earth3DPath ( zenith , azimuth , filename); }
#ifndef EARTH3DMODEL_H
#define EARTH3DMODEL_H

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


//CERN ROOT
#include "TCanvas.h"
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoVolume.h"

class Earth3DModel
{

  public:

  std::string premtable; //Prem .txt file from OscProb/PremTables
  double zenith = TMath::Pi(); //By default it considers a down going neutirno
  double azimuth = 0.0;
  //double Det[3]= {0.0,0.0,-6371.0}; //IceCube site is the Default location
  double Det[3]= {0.0,0.0,-6368.0}; //km3net site is the Default location

  //Specific layers
  int PremRow=1; // Row index in the .txt Prem File
  double LayerDensityContrats=0; // Percentage (%) difference in density  
  double LayerZoAContrats=0; // Percentage (%) difference in zoa (composition) 

  void SetLayerProp(int LayerNumber, double DensityContrats, double ChemicalContrats); //Modify properties of layers

  std::vector< std::vector<double> >  ChangeLayerProp(std::vector< std::vector<double> > EarthMatrix); //Change properties in layer

  //LLVPs
  bool Pile = false; //By Default, model does not include LLVPs
  std::string PileShape = "pancake"; //LLVP shape: "pancake" or "cake"
  double PileThickness = 1000;
  double aWidth = 25.0; //LLVP Angular With
  double PileDensityContrats = 2; // 2% more dense
  double PileZoAContrats = 0.0; // 0% Chemical difference
  double Baseline_LLVP = 0.0; 

  void SetPile( bool value, std::string shape, double density, double zoa ); // Activate the LLVPs
  void CreatePanCake(std::vector<TGeoVolume*> LAYER, std::vector< std::vector<double> > PremMatrix);
  void CreateCake(std::vector<TGeoVolume*> LAYER, std::vector< std::vector<double> > PremMatrix );




  //EARTH MODEL 
  TCanvas *EarthCanvas;

  void SetModel(std::string path2table);

  void SetDetector(double Position[3]);

  void SetDirection(double theta , double phi ); 

  std::vector< std::vector<double> > GetPremData( std::string path2table); // Convert prem .txt file into readable  matrix

  std::vector<std::vector<double>> Earth3DPath ( double zen , double azi, std::string path2table); //Create the 3D Model and calculate Neutrino Path

  std::vector<std::vector<double>> Create3DPath (); // From 3D model, provides the corresponding Neutrino Path inside the Earth

  //std::vector<std::vector<double>> CreateFlickerPath();//To be eliminated


};

#endif
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

class Earth3DModel
{

  public:

  std::string filename; //Prem file
  double zenith;
  double azimuth;
  double Det[3]= {0.0,0.0,-6371.0}; //South Pole is the Default location

  bool Anomaly = false; //By Default we dont need LLVPs
  std::vector<int> WhichLayersLLVPs{1}; //Specify which layers contain an LLVPs segment
  double aWidth = 180.0; //LLVP Angular With
  double drho = 3; // 3% more dense
  double dzoa = 0.0; // 0% Chemical difference
  
  // std::vector< std::vector < double > > TheNuPath;

  void SetModel(std::string model);

  void SetDetector(double Position[3]);

  void SetDirection(double theta , double phi ); 

  
  //LLVPs

  void ActiveHeterogeneity( bool value ); // Activate the LLVPs

  /*

  void SetLLVPAtt( std::vector<int> layers = {10}, double AngularWidth = 45.0, double RhoDiff = 3.0 , double ChemDiff = 0.0) 
  { 

    LLVPIdLayers = layers; //Layers to be modified
    aWidth = AngularWidth;
    drho = 1 + (RhoDiff)/100;
    dzoa = 1 + (ChemDiff)/100;

  }

  */

  std::vector< std::vector<double> > GetPremData( std::string PREM_MODEL = "prem_15layers.txt" );

  int LabelLayer (double radius);

  std::vector<std::vector<double>> Earth3DPath ( double zen , double azi, std::string MODEL);

  std::vector<std::vector<double>> Create3DPath ();


};

#endif
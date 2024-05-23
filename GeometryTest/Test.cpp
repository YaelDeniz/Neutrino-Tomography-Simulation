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


int main()
{

  Earth3DModel test;

  test.SetModel("prem_44layers.txt");

  test.SetDirection(140.0, 0.0);

  //test.LLVPIdLayers()


  std::vector<int> LLVPSegments {25,26,27,28,29,30,31,44};

  test.LLVPIdLayers = LLVPSegments;

  test.aWidth = 70;
  

  //void SetLLVPAtt( LLVPSegments, aWidth, 3.0, 0);
  
  test.ActiveHeterogeneity( true );

  std::vector<std::vector<double>> EarthPath = test.Create3DPath();
//  std::vector<std::vector<double>> EarthPath = CreateEarth3DPath( 180 , 0,  "prem_15layers.txt" );

  double TestL = 0  ;

  //double Ltot = -2.0*6371*cos(th);

  for (int i = 0; i < EarthPath.size(); ++i)
  {


    std::cout << EarthPath[i][0] << " " << EarthPath[i][1] << " " << EarthPath[i][2] << " " << std::endl;
    //TestL = TestL + EarthPath[i][0];
    
  }


return 0;

}
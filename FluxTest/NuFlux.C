#include "NuFlux.h"

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

//ROOT

#include "TH2.h"

using namespace std;


bool NuFlux::DataHeader (const std::string &aline) 
{
    bool aString;

    if (aline.find("flux") != string::npos || aline.find("Enu(GeV)") != string::npos )
    {
        aString = true;
    }

    else {aString = false;}

    return aString;

}

std::vector<std::vector<double>> NuFlux::reshape(const std::vector<double>& vec, int a, int b) {
    // Check if the total number of elements matches
    if (vec.size() != a * b) {
        throw std::invalid_argument("The size of the vector does not match the dimensions of the matrix.");
    }

    // Initialize the reshaped matrix
    std::vector<std::vector<double>> matrix(a, std::vector<double>(b));

    // Fill the matrix with elements from the vector
    for (int i = 0; i < a; ++i) {
        for (int j = 0; j < b; ++j) {
            matrix[i][j] = vec[i * b + j]; //Row Major Order
        }
    }

    return matrix;
}


//-------------------------------------------------------------------------------------------------------------------------------
TH2D* NuFlux::GetFluxHist(int flavor, std::vector<std::vector<double>> FluxData)
{

    std::vector<double> nuflux;

    for (int i = 0; i < FluxData.size(); ++i)
    {
        nuflux.push_back(FluxData[i][flavor]); // {log10(Enu),log10(NuMu*1e-4),log10(NuMubar*1e-4),log10(NuE*1e-4),log10(NuEbar*1e-4)}
    }

    //Create histogram

    //int xbins = sizeof(Log10E) / sizeof(Log10E[0]);
    //int ybins = sizeof(CZ) / sizeof(CZ[0]);

    int xbins = 41;
    int ybins = 20;

    TH2D *flux = new TH2D("flux","flux",xbins,-0.025,2.025,ybins,-1.0,1.0);

    //TGraph2D *flux2D = new TGraph2D(logEpts.size()*CosZpts.size());

    std::vector<std::vector<double>> HistContent =  reshape(nuflux,20,41);

    std::cout << HistContent.size() << " " << HistContent[0].size() << std::endl;

    for (int j = 1; j <=xbins; ++j)
    {
        //int m = j - 1 ;

       // std::cout << "[" <<j <<","<<m <<"]" << std::endl;

        
        for (int i = ybins; i >= 1 ; -- i) //int i = cbins; i >= 1 ; --i
        {
            //std::cout << flux->GetXaxis()->GetBinCenter(j) << " " << flux->GetYaxis()->GetBinCenter(i) << "[" <<j <<","<<m <<"]" << std::endl;
            int m = j - 1 ;
            
            int n =(ybins) - i ;

         //   std::cout << "{" << i << "," << n << "}" << std::endl;

            flux->SetBinContent(j,i,HistContent[n][m]);

            
        }
       // std::cout << " " << std::endl;
    }
 //delete HistContent;
 
 //delete nuflux;

 return flux;
}
//------------------

 std::vector< std::vector<double> > NuFlux::SetFluxData( std::string nuFluxFile = "./SP_AziAveraged_solmin/spl-nu-20-01-000.d" )
{


    std::ifstream HONDA_FLUX;

    //Density profile PREM altered--------------------------------------------------------------------------

    //HONDA_FLUX.open("./AziAveraged-solmin/grn-nu-20-01-000.d"); //km3net

    
    HONDA_FLUX.open(nuFluxFile); //ICE CUBE

        //Parameters for looping the Honda File

    std::string nline;
    
  

    // Check if the file is Open Properly
    
    if ( HONDA_FLUX.is_open() ) { std::cout<< "Flux file opened"<< std::endl; }

    else { std::cout << "Error opening file"<<std::endl; }

    std::vector< std::vector<double> > FluxData;
    
    while(getline(HONDA_FLUX, nline))
    {

        if (!DataHeader(nline)) 
        {   


            // String stream for parsing the input string
            std::istringstream DataRow(nline);
        
            // Vector to store the parsed doubles
            //std::vector<double> DataRow;
        
            // Temporary variable to store each parsed double
            double  Enu,NuMu,NuMubar,NuE,NuEbar;
        
            // Parse the input string and store the doubles in the vector
            DataRow >> Enu >> NuMu >> NuMubar >> NuE >> NuEbar;


            if (Enu>= 1.0 && Enu<= 100.0 )
            {
              
         
                FluxData.push_back( {log10(Enu),log10(NuMu*1e-4),log10(NuMubar*1e-4),log10(NuE*1e-4),log10(NuEbar*1e-4)} );

             //   std::cout << Enu <<" "<< NuMu << " " << NuMubar << " " << NuE << " " << NuEbar << std::endl;

            }

        }

        else 
        { 
            std::cout << " " << std::endl;
         }      
        
    }

    return FluxData;
}





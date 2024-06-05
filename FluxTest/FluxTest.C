
bool DataHeader (const std::string &aline) 
{
    bool aString;

    if (aline.find("flux") != string::npos || aline.find("Enu(GeV)") != string::npos )
    {
        aString = true;
    }

    else {aString = false;}

    return aString;

}

std::vector<std::vector<double>> reshape(const std::vector<double>& vec, int a, int b) {
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




void FluxTest()
{


    std::ifstream HONDA_FLUX;

    //Density profile PREM altered--------------------------------------------------------------------------

    HONDA_FLUX.open("./AziAveraged-solmin/grn-nu-20-01-000.d");

        //Parameters for looping the Honda File

    std::string nline;
    
  

    // Check if the file is Open Properly
    
    if ( HONDA_FLUX.is_open() ) { std::cout<< "Flux file opened"<< std::endl; }

    else { std::cout << "Error opening file"<<std::endl; }

    double  Enu,NuMu,NuMubar,NuE,NuEbar;

    std::vector< std::vector<double> > FluxData;

   

    bool AXIS = true;
    std::vector < double > l10E; 

    double CosZ = -1.0;
    while(getline(HONDA_FLUX, nline))
    {

        if (!DataHeader(nline)) 
        {   


            // String stream for parsing the input string
            std::istringstream iss(nline);
        
            // Vector to store the parsed doubles
            std::vector<double> DataRow;
        
            // Temporary variable to store each parsed double
            double  Enu,NuMu,NuMubar,NuE,NuEbar;
        
            // Parse the input string and store the doubles in the vector
            iss >> Enu >> NuMu >> NuMubar >> NuE >> NuEbar;


            if (Enu>= 1.0 && Enu<= 100.0 )
            {
              
                //DataRow.push_back(Enu);
                //DataRow.push_back(NuMu);
                //DataRow.push_back(NuMubar);
                //DataRow.push_back(NuE);
                //DataRow.push_back(NuEbar);
                FluxData.push_back( {CosZ,log10(Enu),log10(NuMu*1e-4),log10(NuMubar*1e-4),log10(NuE*1e-4),log10(NuEbar*1e-4)} );

             //   std::cout << Enu <<" "<< NuMu << " " << NuMubar << " " << NuE << " " << NuEbar << std::endl;

            }

        }

        else { CosZ += 0.05;
         std::cout << " " << std::endl; }      
        
    }


    //Sort Angular Bins
    std::vector<double> CosZpts;

    for (int i = 0; i < FluxData.size(); ++i)
    {
        CosZpts.push_back(FluxData[i][0]-0.1/2.0);
    }

    //CosZpts.push_back(1.1);

    set<double> s1;
    unsigned size1 = CosZpts.size();
    for( unsigned i = 0; i < size1; ++i ) s1.insert( CosZpts[i] );
    CosZpts.assign( s1.begin(), s1.end() );
    double CZ[CosZpts.size()];
    copy(CosZpts.begin(),CosZpts.end(),CZ);


//    for (int i = 0; i < CosZpts.size(); ++i)
  //  {
     //   std::cout << CZ[i]<<std::endl;
   // }





    //Sort Energy Bins
    
    std::vector<double> logEpts;

    for (int i = 0; i < FluxData.size(); ++i)
    {
        logEpts.push_back(FluxData[i][1]);
    }

    set<double> s2;
    unsigned size = logEpts.size();
    for( unsigned i = 0; i < size; ++i ) s2.insert( logEpts[i] );
    logEpts.assign( s2.begin(), s2.end() );
    double Log10E[logEpts.size()];
    copy(logEpts.begin(),logEpts.end(),Log10E);


/*
    for (int i = 0; i < logEpts.size(); ++i)
    {
        std::cout << " * " <<Log10E[i]<<std::endl;
    }

*/
    std::vector<double> nuflux;
    for (int i = 0; i < FluxData.size(); ++i)
    {
        nuflux.push_back(FluxData[i][4]); // {CosZ,log10(Enu),log10(NuMu*1e-4),log10(NuMubar*1e-4),log10(NuE*1e-4),log10(NuEbar*1e-4)}
    }

    //Matrix for Histogram & Histogram Draw


    //Create histogram

    //int xbins = sizeof(Log10E) / sizeof(Log10E[0]);
    //int ybins = sizeof(CZ) / sizeof(CZ[0]);

    int xbins = 41;
    int ybins = 20;


    TCanvas *c1 = new TCanvas();
    TH2D *flux = new TH2D("flux","flux",xbins,-0.025,2.025,ybins,-1.0,1.0);

    //TGraph2D *flux2D = new TGraph2D(logEpts.size()*CosZpts.size());

    std::vector<std::vector<double>> HistContent =  reshape(nuflux,CosZpts.size(), logEpts.size());

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

            std::cout << "{" << i << "," << n << "}" << std::endl;

            flux->SetBinContent(j,i,HistContent[n][m]);

            
        }
        std::cout << " " << std::endl;
    }

    flux->SetStats(0);

flux-> SetTitle("FLux");

flux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
flux->GetYaxis()->SetTitle("Cos(#theta_{z})");
flux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");


flux->Draw("COLZ");







    double Emax = 100;

    double Emin = 1;

    double pts = 1000; 

    double hE = (100 - 1)/pts;

    double Ei = 1;

    TCanvas *c2 = new TCanvas();

    TGraph  *test = new TGraph(pts);

 std::cout << "----------------" << std::endl;

    for (int i = 0; i < pts; ++i)
     {
          

         //double logphi = flux->Interpolate(log10(Ei),-0.9);

         //std::cout << Ei + i*hE << " " << Ei << " " << log10(Ei + i*hE) << " " << pow(10, log10(Ei + i*hE)) << std::endl;

         double x =  Ei + i*hE;

         double logx = log10(x);
         
         double logphi = flux->Interpolate(logx,-0.9);

         //std::cout <<  logphi << std::endl;

         double y = pow(10,logphi)*x*x*x;

         test -> SetPoint (i,  logx,  y);
     } 



    

    test->GetXaxis()->SetTitle("log_{10}(E/GeV)");
    
    test->GetYaxis()->SetTitle("#Phi_{#nu}*E_{#nu}^{3}");

    test->SetTitle("");
    
    test->Draw();







    /*

    std::vector<std::vector<double>> HistContent =  reshape(nuflux,CosZpts.size(), logEpts.size());

    //std::cout << HistContent.size() << " " << HistContent[0].size() << std::endl;
    /*
    for (int i = 0; i < 5 HistContent.size() ; ++i)
    {
        std::cout << i <<" :  " ;
        for (int j = 0; j < HistContent[0].size() ; ++j)
        {
            std::cout << HistContent[i][j] << " " ;

            flux->SetBinContent(i,j,HistContent[i][j]);
        }

        std::cout << std::endl;
    }
    

    int k = 0;

    int ebins = HistContent[0].size();

    int cbins = HistContent.size();

    //int ebins = 3;

    //int cbins = 14;

    for (int j = 0; j < ebins ; ++j) //41 >
    {

        int m = (ebins-1) - j;

        std::cout << " " << std::endl;
 
        for (int i = cbins; i >= 1 ; --i) //20 ^
        {
         
            int n =(cbins) - i ;

            std::cout << "Index check " << i << " " << n << std::endl;

            flux->SetBinContent(j,i,HistContent[n][j]);


            //flux2D->SetPoint(k,1,1, 1 );
            k++;
        
        }

    } 

flux->SetStats(0);

flux-> SetTitle("FLux");

flux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
flux->GetYaxis()->SetTitle("Cos(#theta_{z})");
flux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");


flux->Draw("COLZ");

TH1D *Bins = new TH1D("bins","bins",41,-0.025, 2.025 );



for (int i = 1; i <= 41; ++i)
{
    std::cout <<   Bins->GetBinCenter(i) << std::endl; 
    std::cout <<   Bins->GetBinWidth(i) << std::endl; 
    std::cout <<   Bins->GetBinLowEdge(i)<< std::endl;
    std::cout << " " << std::endl;
}


 /*
TCanvas *c3 = new TCanvas();
flux2D->Draw("surf1");


    double Emax = 100;

    double Emin = 1;

    double pts = 1000; 

    double hE = (100 - 1)/pts;

    double Ei = 1;

    TCanvas *c2 = new TCanvas();

    TGraph  *test = new TGraph(pts);

 std::cout << "----------------" << std::endl;

    for (int i = 0; i < pts; ++i)
     {
          

         //double logphi = flux->Interpolate(log10(Ei),-0.9);

         //std::cout << Ei + i*hE << " " << Ei << " " << log10(Ei + i*hE) << " " << pow(10, log10(Ei + i*hE)) << std::endl;

         double x =  Ei + i*hE;

         double logx = log10(x);
         
         double logphi = flux->Interpolate(logx,-0.9);

         //std::cout <<  logphi << std::endl;

         double y = pow(10,logphi)*x*x*x;

         test -> SetPoint (i,  logx,  y);
     } 



    

    test->GetXaxis()->SetTitle("log_{10}(E/GeV)");
    
    test->GetYaxis()->SetTitle("#Phi_{#nu}*E_{#nu}^{3}");

    test->SetTitle("");
    
    test->Draw();




    //Populate histogram

    //std::cout <<  "Matrix Visualization" << " "  << FluxData.size() << " " << FluxData[0].size()   << std::endl;


    //for (int j = 0; j < ybins; ++j)
    //{

      //  flux->SetBinContent(j,j,j) ;
    //}


/*
    for (int j = 0; j < ybins; ++j)
    {
        for (int i = 0; i < xbins; ++i)
        {
            //std::cout <<  FluxData[i][j]<< " " ; 

            std::cout << j << " " << i << " "  << FluxData[i][0] << " "  << FluxData[i][1] << std::endl;


        }
        
        std::cout << std::endl ;  
    }


 */



    




    /*

    set<double> s;
    unsigned size = l10E.size();
    for( unsigned i = 0; i < size; ++i ) s.insert( l10E[i] );
    l10E.assign( s.begin(), s.end() );
    double Log10E[l10E.size()];
    copy(l10E.begin(),l10E.end(),Log10E);


    
    for (int i = 0; i < l10E.size(); ++i)
    {
        std::cout <<"*" << Log10E[i] << " " <<std::endl;
    }

    std::vector < double > zen;

    for (int i = 0; i <= 20; ++i)
    {
        zen.push_back(-1.0+i*0.1);
        std::cout <<  zen[i] << std::endl;
    }

    int xbins = Log10E.size();
    int ybins = zen.size();

    TCanvas *c1 = new TCanvas();
    TH2D *flux = new TH2D("flux","flux",xbins,0.0,2.0,ybins,-1.0,1.0);
  */

}

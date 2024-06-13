class NuFlux
{

    public: 

    //Neutrino Energies
    double Enu_min= 1.0; \\GeV
    double Enu_max= 100.0;

    //Detector Locations
    std::string  DetLoc= "South Pole";

    bool DataHeader (const std::string &aline); // Get data into tabular form but ignore Headers from file

    std::vector<std::vector<double>> reshape(const std::vector<double>& vec, int a, int b); //Reshape colum vector into matrix

    std::vector< std::vector<double> > SetFluxData( std::string nuFluxFile = "./SP_AziAveraged_solmin/spl-nu-20-01-000.d" );
    
    TH2D*GetFluxHist(int flavor, std::vector<std::vector<double>> FluxData); // Create 2D Histogram

    
}


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



//------------------
void NuFlux_ClassTest()
{
 
NuFlux SPflux;

std::vector< std::vector<double> > FluxData = SPflux.SetFluxData();


//Matrix for Histogram & Histogram Draw

//-------------------------------------------------------------------------------------------------------------------------------
TH2D* muflux =  GetFluxHist(1,FluxData); //MuFlux
TH2D* mubflux =  GetFluxHist(2,FluxData); //MuBarFlux
TH2D* eflux =  GetFluxHist(3,FluxData); //EFlux
TH2D* ebflux =  GetFluxHist(4,FluxData); //EBarFlux

//---------------------------------------------------------------------------------------------------------------------------------------
muflux->SetStats(0);
muflux-> SetTitle("#nu_{#mu} Flux");
muflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
muflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//muflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

mubflux->SetStats(0);
mubflux-> SetTitle("#bar{#nu}_{#mu} Flux");
mubflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
mubflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//mubflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

eflux->SetStats(0);
eflux-> SetTitle("#nu_{e} Flux");
eflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
eflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//eflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

ebflux->SetStats(0);
ebflux-> SetTitle("#bar{#nu}_{e} Flux");
ebflux->GetXaxis()->SetTitle("log_{10}(E/GeV)");
ebflux->GetYaxis()->SetTitle("Cos(#theta_{z})");
//ebflux->GetZaxis()->SetTitle("log_{10}(#Phi_{#nu})");

TCanvas *c = new TCanvas();
c->Divide(2,2);

c->cd(1);
eflux->Draw("COLZ");


c->cd(2);
ebflux->Draw("COLZ");

c->cd(3);
muflux->Draw("COLZ");

c->cd(4);
mubflux->Draw("COLZ");

//c->SetCanvasSize(500, 500);

c->SetWindowSize(1200, 800);

//gStyle->SetPadLeftMargin(10000); 
//gStyle->SetPadRightMargin(0.1);




//Interpolation of the flux





    double Emax = 100;

    double Emin = 1;

    double pts = 1000; 

    double hE = (100 - 1)/pts;

    double Eo = 1;

    

    TGraph  *numu = new TGraph(pts);
    TGraph  *numub = new TGraph(pts);
    TGraph  *nue = new TGraph(pts);
    TGraph  *nueb = new TGraph(pts);

 std::cout << "----------------" << std::endl;

    for (int i = 0; i < pts; ++i)
     {
          

 

         double Ei =  Eo + i*hE;

         double logEi = log10(Ei);
         
         double logphimu = muflux->Interpolate(logEi,-0.9);
         double logphimub = mubflux->Interpolate(logEi,-0.9);
         double logphie = eflux->Interpolate(logEi,-0.9);
         double logphieb = ebflux->Interpolate(logEi,-0.9);

         

         double phimu = pow(10,logphimu)*Ei*Ei*Ei;
         double phimub = pow(10,logphimub)*Ei*Ei*Ei;
         double phie = pow(10,logphie)*Ei*Ei*Ei;
         double phieb = pow(10,logphieb)*Ei*Ei*Ei;

         numu -> SetPoint (i,  logEi,  phimu);
         numub -> SetPoint (i,  logEi,  phimub);
         nue -> SetPoint (i,  logEi,  phie);
         nueb -> SetPoint (i,  logEi,  phieb);
     } 



    TCanvas *c1 = new TCanvas();
    TMultiGraph * multi = new TMultiGraph();

    numu->SetLineColor(4);
    nue->SetLineColor(2);
    numub->SetLineColor(4);
    nueb->SetLineColor(2);

    numub->SetLineStyle(9);
    nueb->SetLineStyle(9);

    multi->Add(numu);
    multi->Add(numub);
    multi->Add(nue);
    multi->Add(nueb);



    multi->GetXaxis()->SetLimits(0.0,2.0);
    multi->GetXaxis()->SetTitle("log_{10}(E/GeV)");

    multi->GetYaxis()->SetLimits(0.0,0.03);
    multi->GetYaxis()->SetTitle("#Phi_{#nu}*E^{3}");

    multi->SetTitle("Bilinear interpolation of Honda flux");

    gPad->SetGrid(1,1);
    multi->Draw("A");

    c->SaveAs("./FluxHistSP.jpg");
    c1->SaveAs("./FluxInterpolSP.jpg");

   
   
}

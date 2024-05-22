std::vector< std::vector<double> > GetPremData( std::string PREM_MODEL = "prem_15layers.txt" )
{
   std::string PREM_PATH = "/home/dehy0499/OscProb/PremTables/"+PREM_MODEL;

   //--- open PREM model 
   std::vector< std::vector<double> > PremMatrix; // Matrix data  form of Prem tables
   std::vector<double> PremRow(4);               // [radius density z/a, layer ID]

   // Variables for storing table rows
   double radius, density, zoa, layer;

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

int LabelLayer (double radius)
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







void BinnedEarth()
{
  //Read Prem Data

  gSystem->Load("libGeom");

  std::vector< std::vector<double> > PremMatrix = GetPremData(); // Sort PREM model into a readable matrix

  std::vector< std::vector<double> > LLVPMatrix = GetPremData(); // A Copy of PREM model matrix to construct LLVPs

  



  double rPREM = PremMatrix.back()[0]; // Outermost Layer in the PREM model

  std::vector< std::vector<double> > EarthPath;  // Store Paths inside the Earth for oscillation calculations

  std::vector< int > EarthPathId;  // Store ID layer for each Paths inside the Earth for oscillation calculations



  

  TCanvas *c1 = new TCanvas("3D Earth", "3D Earth",0,0,600,600);
 
  new TGeoManager("EarthModel3D", "Simple 3D geometry");

  //Define world (Volume container)

  //Defining Vacuum Medium
  TGeoMaterial *VAC = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
  
  TGeoMedium *vac = new TGeoMedium("VACUUM",1,VAC);
  
  TGeoVolume *top = gGeoManager->MakeBox("TOP",vac,rPREM,rPREM,rPREM);

  gGeoManager->SetTopVolume(top);

  //Create Earth

  /* Each Layer is represented of a Spherical Shell defined by a Inner Radius (rmin) and Outer Radius(rmax).
  Each Layer Medium is defined by a Material, a Material is created from a specific density value and Z[proton]/A[neutron + proton] (zoa)*/
  double rmin, rmax, density, zoa, layer; 

  std::vector<TGeoMaterial*> MAT; //Vector to define material for each layer en in PREM (Density & Z/A)
  std::vector<TGeoMedium*>  MED;  //Vector to define medium for each layer en in PREM 
  std::vector<TGeoVolume*> LAYER; //Vector to define a  volume for each layer en in PREM

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

      LAYER.push_back( gGeoManager->MakeSphere(LayerName,MED[i], rmin,rmax,0,180,0,360) ); // Define Volume for the ith Layer

      LAYER[i]->SetVisibility(kFALSE);
  }
 

  //Definition LLVPS

  /*In this version of the Code LLVPs are just segmets of SPHERE inside specific layers in the Lower Mantle*/

  std::vector<TGeoMaterial*> LLVPMat; //Vector to define material for each LLVP segment ( A X% anomaly in local Density or Z/A)
  std::vector<TGeoMedium*> LLVPMed;   //Vector to define medium for each LLVP segment
  std::vector<TGeoVolume*> LLVPLAYER; //Vector to define Volume for each LLVP segment
  
  TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);// Some Geometrical trasformation that move LLVPs to the right place

  int LLVPIdLayers[1] = {7}; //Layers to be modified

  int LLVPint = sizeof(LLVPIdLayers) / sizeof(int);

  double gamma = 45.0/LLVPint; // Angular width of each LLVP segment
  
  double drho, dzoa, gammai, rminLLVP, rmaxLLVP;
 
  for (int i = 0; i < LLVPint; ++i)
  {

    int index = LLVPIdLayers[i]-1; //Index of the layer that will contain LLVPs

    //Generate segment names
    char integer_string[32];
    
    sprintf(integer_string, "%d", i+1);

    char llvpmaterial_string[64]="LLVPMaterial";
    strcat(llvpmaterial_string, integer_string);
    const char *llvpMatName = llvpmaterial_string;

    char llvpmedium_string[64]="LLVPMedium";
    strcat(llvpmedium_string, integer_string);
    const char *llvpMedName = llvpmedium_string;
    
    char llvplayer_string[64]="LLVPLayer";
    strcat(llvplayer_string, integer_string);
    const char *llvpLayerName= llvplayer_string;

    rminLLVP = LLVPMatrix[index-1][0];
    rmaxLLVP = LLVPMatrix[index][0];
    drho = (1.03)*LLVPMatrix[index][1];
    dzoa = (1.00)*LLVPMatrix[index][2];
    gammai =(LLVPint - i)*gamma;
    LLVPMat.push_back( new TGeoMaterial(llvpMatName, 1, dzoa, drho) ); // Create LLVP material
  
    LLVPMed.push_back( new TGeoMedium(llvpMedName,1,LLVPMat[i]) ); // Create LLVP medium

    LLVPLAYER.push_back( gGeoManager -> MakeSphere(llvpLayerName, LLVPMed[i] , rminLLVP , rmaxLLVP , 0,gammai,0,360 ) ); // DEFINE LLVP Segment

    LLVPLAYER[i]->SetLineColor(kGreen);

    LAYER[index]->AddNode(LLVPLAYER[i],1,rot1);

  }
  
  //Add Earth Layers to TOP Volume

  for (int i = 0; i < PremMatrix.size(); ++i)
  {
        top->AddNode(LAYER[i],1);
  }

  gGeoManager->CloseGeometry(); // Finish Geometry

  /*
  
  gGeoManager->SetTopVisible(); // the TOP is invisible
 
  top->Draw();

  TView *view = gPad->GetView();
  
  view->ShowAxis();

  */
 
  // Tracking----------------------------------------------------------------------------------------------------------------------------------------

  //GetInitialPosition(zen,azi)
  //Calculate Paths inside the Earth

   //Direction of neutrino in spherical coordiates: https://mathworld.wolfram.com/SphericalCoordinates.html
  double zen = 180.0; // zenith angle respect to the detector location in degrees (90,180]

  double azi = 0.0;  //Azimuthal angle respect to the detector location in degrees [0,360]

  double th = TMath::Pi()*( 1-(zen/180.0) ); //Polar angle in spherical coordinates in rads [pi/2, pi] 

  double phi = TMath::Pi()*(azi)/180.0; //Azimuthal angle in spehericla coordinates in rads [0, 2*pi]
   
   
  double r = rPREM; // Max Radius

  double Det[3]= {0.0,0.0,-6371.0}; //Detector location in cartesian coordinates
  
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

  //TPolyMarker3D *l1 = new TPolyMarker3D(2,2); //Markers indicating DetLoc and Inical Neutrino
  //l1->SetPoint( 0 , o[0],o[1],o[2]); //Detector
  //l1->SetPoint( 1 ,xo,yo,zo); //Incoming Neutrino

  //TPolyLine3D *l2 = new TPolyLine3D(); // Lines that represent neutrino Paths.
  


  int i = 0;

  double Dx2, Dy2, Dz2, Dnorm;

  double sumL = 0; //Track Total Baseile

  //Neutrino Propagation inside the Earth
  while (!gGeoManager->IsOutside ())
  {  

         std::cout << "Starting loop" << std::endl;
         
         //gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
      
         int nodeID= gGeoManager->GetCurrentNode()->GetIndex();

         const Double_t *cpoint = gGeoManager->GetCurrentPoint(); // Current nuetrino Position

         const char *path = gGeoManager->GetPath();

         TGeoVolume *cvol = gGeoManager->GetCurrentVolume();

         TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial(); //Material of the current Boundary
            

         //Double_t Li = gGeoManager->GetStep(); //Baseline Segement

         Double_t safety = gGeoManager->GetSafeDistance(); //Distance to the next boundary/Neutrino Baseline

         gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
          
         Double_t Li = gGeoManager->GetStep(); //Baseline Segement

         double R_i = sqrt(cpoint[0]*cpoint[0]+ cpoint[1]*cpoint[1] + cpoint[2]*cpoint[2]); //Current Radius

         EarthPath.push_back({Li, cmat->GetDensity(), cmat->GetZ(),R_i}); //Store Path information [Baseline segment, Density, Z/A, R_i]
        
         EarthPathId.push_back(LabelLayer(R_i)); // Store layer ID

         sumL = sumL + Li; //Track Total Segment (Max = -2*R_earth*cos(zen))

         // l2->SetPoint( i , cpoint[0],cpoint[1],cpoint[2]);




        //Path Infromation 

         std::cout << "Current path is: " << path << " Layer radius: "  << R_i << " id: " << LabelLayer(R_i) <<  std::endl;

         std::cout << "L_i " << Li  << " rho "  << cmat->GetDensity() << " zoa " << cmat->GetZ() << " medId " << cvol->GetMedium()->GetId() << std::endl;

         std::cout << std::endl; 

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

  //l1->SetMarkerColor(6);
  //l1->SetMarkerSize(3);
  //l2->SetLineColor(14);
  //l2->SetLineWidth(3);
   
  //l1->Draw("same");
  //l2->Draw("same");

  //DrawPrem

  //TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);




  for (int i = 0; i < EarthPath.size(); ++i)
  {
    std::cout << EarthPath[i][0] << " " << EarthPath[i][1] << " " << EarthPath[i][2] << " " << EarthPath[i][3] << " " << EarthPathId[i] << std::endl;
  }

/*
// Number of layers
  int nlayers = EarthPath.size();

  // Create a TGraph to plot
  TGraph* gr = new TGraph(nlayers);
  TGraph* prem = new TGraph(PremMatrix.size());

  // Loop over paths
  for(int i=1; i<16; i++){
    gr->SetPoint(i, EarthPath[i][3], EarthPath[i][1]);

    std::cout << EarthPath[i][3] << " " << EarthPath[i][1] << std::endl; 
   
  }

  for(int i=0; i < PremMatrix.size(); i++){
   
    prem->SetPoint((PremMatrix.size()-1)-i, PremMatrix[i][0], PremMatrix[i][1]);

    std::cout << " * " <<PremMatrix[i][0] << " " << PremMatrix[i][1] << std::endl; 
  }
 
  gr->SetLineColor(kRed);
  
  gr->Draw("C*");


  prem->SetLineWidth(3);
  prem->SetMarkerStyle(21);
  prem->SetLineColor(6);
  prem->Draw("C");
  */ 





  ////



   //Geometrical display
   std::cout << "Detector Geometrical settings  d1 " << d1 << " d2 " << d2 << std::endl;
   std::cout << " Detector: "  << Det[0] << " " << Det[1] << " " << Det[2] << std::endl;
   std::cout << " Neutrino: "  << xo << " " << yo << " " << zo << std::endl;
   std::cout << " Neutrino: "  << sqrt(xo*xo +yo*yo + zo*zo) << std::endl;

}
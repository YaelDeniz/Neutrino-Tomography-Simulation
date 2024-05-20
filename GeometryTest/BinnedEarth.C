std::vector< std::vector<double> > GetPremData( std::string PREM_MODEL = "prem_44layers.txt" )
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

  std::vector< std::vector<double> > PathMatrix;
 
  
  std::vector<double> PathData(4); // Row [L, rho, z/a, layer]

  std::vector< std::vector<double> > PremMatrix = GetPremData(); // Sort PREM model into a readable matrix

  std::vector< std::vector<double> > LLVPMatrix = GetPremData();

  



  double rEarth = PremMatrix.back()[0];

  std::vector< std::vector<double> > EarthPath;
  std::vector< int > EarthPathId;



  

  TCanvas *c1 = new TCanvas("c", "c",0,0,600,600);
 
  new TGeoManager("simple1", "Simple geometry");

  //Define world
  TGeoMaterial *VAC = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
  
  TGeoMedium *vac = new TGeoMedium("VACUUM",1,VAC);
  
  TGeoVolume *top = gGeoManager->MakeBox("TOP",vac,rEarth+10,rEarth+10,rEarth+10);

  gGeoManager->SetTopVolume(top);

  //Create Earth
  std::vector<TGeoMaterial*> MAT; //
  std::vector<TGeoMedium*>  MED; //
  std::vector<TGeoVolume*> LAYER; //

  double rmin, rmax, density, zoa, layer;

  for (int i = 0; i < PremMatrix.size(); ++i)
   {
      
      if (i==0)
      {
        rmin = 0;
      }

      else
      {

      rmin = PremMatrix[i-1][0];
      
      }

      rmax = PremMatrix[i][0];

      density = PremMatrix[i][1]; 
      
      zoa = PremMatrix[i][2];
      
      layer = PremMatrix[i][3];

      //Generate Names
      char integer_string[32];
      sprintf(integer_string, "%d", i+1);

      char material_string[64]="Material";
      char medium_string[64]="Medium";
      char layer_string[64]="Layer";

      strcat(material_string, integer_string);
      strcat(medium_string, integer_string);
      strcat(layer_string, integer_string);
      

      const char *MatName = material_string;
      const char *MedName = medium_string;
      const char *LayerName= layer_string;

      

      MAT.push_back( new TGeoMaterial(MatName,1.0 ,zoa,density) );
       
      MED.push_back( new TGeoMedium(MedName,1, MAT[i]) );

      LAYER.push_back( gGeoManager->MakeSphere(LayerName,MED[i], rmin,rmax,0,180,0,360) );//0
      LAYER[i]->SetVisibility(kFALSE);

  
      //top->AddNode(LAYER[i],1);
  }

  LAYER[23]->SetVisibility(kTRUE);
  //LAYER[35]->SetVisibility(kTRUE);
  //LAYER[41]->SetVisibility(kTRUE);



  //Add LLVPS

  std::vector<TGeoMaterial*> LLVPMat;
  std::vector<TGeoMedium*> LLVPMed;
  std::vector<TGeoVolume*> LLVPLAYER;
  
  TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);//Rotation

  int LLVPIdLayers[3] = {25,26,27}; //Layers to be modified

  int LLVPint = sizeof(LLVPIdLayers) / sizeof(int);

  double angle = 45.0/LLVPint;

  double per = 1.03; //LLVP 3% more dense
 
  for (int i = 0; i < LLVPint; ++i)
  {

    int index = LLVPIdLayers[i]-1;

    //Generate Names
    char integer_string[32];
    sprintf(integer_string, "%d", i+1);

    char llvpmaterial_string[64]="LLVPMaterial";
    char llvpmedium_string[64]="LLVPMedium";
    char llvplayer_string[64]="LLVPLayer";

    strcat(llvpmaterial_string, integer_string);
    strcat(llvpmedium_string, integer_string);
    strcat(llvplayer_string, integer_string);
    

    const char *llvpMatName = llvpmaterial_string;
    const char *llvpMedName = llvpmedium_string;
    const char *llvpLayerName= llvplayer_string;
     
    LLVPMat.push_back( new TGeoMaterial(llvpMatName, 1, LLVPMatrix[index][2], per*LLVPMatrix[index][1]) );
  
    LLVPMed.push_back( new TGeoMedium(llvpMedName,1,LLVPMat[i]) );

    LLVPLAYER.push_back( gGeoManager -> MakeSphere(llvpLayerName, LLVPMed[i] ,LLVPMatrix[index-1][0],LLVPMatrix[index][0] , 0,(3-i)*angle,0,360 ) );

    LLVPLAYER[i]->SetLineColor(kGreen);

    LAYER[index]->AddNode(LLVPLAYER[i],1,rot1);
  }
  
  //Add Earth Layers

  for (int i = 0; i < PremMatrix.size(); ++i)
  {
        top->AddNode(LAYER[i],1);
  }

  gGeoManager->CloseGeometry();
  gGeoManager->SetTopVisible(); // the TOP is invisible

  top->Draw();
  TView *view = gPad->GetView();
  view->ShowAxis();
 
  // Tracking----------------------------------------------------------------------------------------------------------------------------------------

  //GetInitialPosition(zen,azi)
  //Calculate Paths inside the Earth

   //Direction of neutrino in spherical coordiates: https://mathworld.wolfram.com/SphericalCoordinates.html
  double zen = 145.0; // zen (90,180]

  double azi = 0.0;

  double th = TMath::Pi()*( 1-(zen/180.0) ); // 0 to Pi

  double phi = TMath::Pi()*(azi)/180.0; // 0 to 2*Pi
   
  //Initial Position/ Direction Cosines
  double xn = sin(th)*cos(phi);
  double yn = sin(th)*sin(phi);
  double zn = cos(th);
   
  double r = rEarth;

  //double DETR = 6368;
  //double r = DETR;

  ROOT::Math::SVector<double, 3> u(xn,yn,zn); // Direction of neutrino
  ROOT::Math::SVector<double, 3> o(0.0,0.0,-r); // Origin of the line
  ROOT::Math::SVector<double, 3> cc(0,0,0); // Origin of the Sphere 
 

  double a = ROOT::Math::Dot( u, u);
  double b = 2*ROOT::Math::Dot( u, o - cc);
  double c = ROOT::Math::Dot( o-cc, o-cc) - r*r ;
  double d1 = (-b+sqrt(b*b -4*a*c))/(2*a);
  double d2 = (-b-sqrt(b*b -4*a*c))/(2*a);

  std::cout << " d1 " << d1 << " d2 " << d2 << std::endl;

  TPolyLine3D *l1 = new TPolyLine3D();
  TPolyLine3D *l2 = new TPolyLine3D();
   //l->SetPoint( 1 , o[0]-0.5*d1*u[0],o[0]-0.5*d1*u[0],o[0]-0.5*d1*u[0]);

  // Initial Neutrino Positionin the Earth

  double xo = o[0]+d1*u[0];
  double yo = o[1]+d1*u[1];
  double zo = o[2]+d1*u[2];

  // Neutrino Direction

  double ni = -u[0];
  double nj = -u[1];
  double nk = -u[2];
  
  gGeoManager->InitTrack(xo, yo, zo, ni, nj, nk); // d*u is the point in the sphere
  
  l1->SetPoint( 0 , o[0],o[1],o[2]);
  
  l1->SetPoint( 1 ,xo,yo,zo);

  int i = 0;

  while (!gGeoManager->IsOutside ())
  {  

         //gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
      
         int nodeID= gGeoManager->GetCurrentNode()->GetIndex();

         const Double_t *cpoint = gGeoManager->GetCurrentPoint();

         const char *path = gGeoManager->GetPath();

         TGeoVolume *cvol = gGeoManager->GetCurrentVolume();

         TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial(); //Material of the current Boundary
            

         Double_t Li = gGeoManager->GetStep(); //Baseline Segement

         Double_t safety = gGeoManager->GetSafeDistance(); //Distance to the next boundary/Neutrino Baseline

         double R_i = sqrt(cpoint[0]*cpoint[0]+ cpoint[1]*cpoint[1] + cpoint[2]*cpoint[2]);
         

         std::cout << "Current path is: " << path << " r "  << R_i << "id: " << LabelLayer(R_i) <<  std::endl;

         std::cout << "L " << Li  << " rho "  << cmat->GetDensity() << " zoa " << cmat->GetZ() << " "<< cmat->GetA() << " " << cvol->GetMedium()->GetId() << std::endl;

         std::cout << std::endl; 
         
         l2->SetPoint( i , cpoint[0],cpoint[1],cpoint[2]);

         gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
         
         i += 1;

        
         EarthPath.push_back({Li, cmat->GetDensity(), cmat->GetZ()});
         EarthPathId.push_back(LabelLayer(R_i));
  }


  for (int i = 0; i < EarthPath.size(); ++i)
  {
    std::cout << EarthPath[i][0] << " " << EarthPath[i][1] << " " << EarthPath[i][2] << " " << EarthPathId[i] << std::endl;
  }

   l1->SetLineColor(6);
   l1->SetLineWidth(3);
   l2->SetLineColor(14);
   l2->SetLineWidth(3);
   
   l1->Draw("same");
   l2->Draw("same");


}
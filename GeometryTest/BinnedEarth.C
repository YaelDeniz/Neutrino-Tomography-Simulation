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




void BinnedEarth()
{
  //Read Prem Data

  gSystem->Load("libGeom");

  std::vector< std::vector<double> > PathMatrix;
  
  std::vector<double> PathData(4); // Row [L, rho, z/a, layer]

  std::vector< std::vector<double> > PremMatrix = GetPremData(); // Sort PREM model into a readable matrix

  double rEarth = PremMatrix.back()[0];

  

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

  
      top->AddNode(LAYER[i],1);
  }

  gGeoManager->CloseGeometry();
  gGeoManager->SetTopVisible(); // the TOP is invisible

  top->Draw();
  TView *view = gPad->GetView();
  view->ShowAxis();
 
  // Tracking----------------------------------------------------------------------------------------------------------------------------------------

  //Calculate Paths inside the Earth

   //Direction of neutrino in spherical coordiates: https://mathworld.wolfram.com/SphericalCoordinates.html
  double zen = 180.0; // zen (90,180]

  double th = TMath::Pi()*( 1-(zen/180.0) ); // 0 to Pi
   
  double phi = 0.0; // 0 to 2*Pi
   
   //Initial Position
  double xn = sin(th)*cos(phi);
  double yn = sin(th)*sin(phi);
  double zn = cos(th);

   //double n = sqrt(x*x + y*y + z*z);

   //Unit vector
   //double nx = x/n;
   //double ny = y/n;
   //double nz = z/n;

   //Origin of the Sphere

   //double zo = Ro/Ro;

   //Initial Point must be in the outermost layer -> Intersection Line spehere

   
  double r = rEarth;

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
  l1->SetPoint( 0 , 0.0, 0.0, -r);
  l1->SetPoint( 1 , o[0]+d1*u[0],o[1]+d1*u[1],o[2]+d1*u[2]);
   //l->SetPoint( 1 , o[0]-0.5*d1*u[0],o[0]-0.5*d1*u[0],o[0]-0.5*d1*u[0]);
  
   
   

  gGeoManager->InitTrack(o[0]+d1*u[0], o[1]+d1*u[1], o[2]+d1*u[2], -u[0] , -u[1], -u[2]); // d*u is the point in the sphere

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
         

         std::cout << "Current path is: " << path << " r "  << sqrt(cpoint[0]*cpoint[0]+ cpoint[1]*cpoint[1] + cpoint[2]*cpoint[2])<< std::endl;

         std::cout << "L " << Li  << " rho "  << cmat->GetDensity() << " zoa " << cmat->GetZ() << " "<< cmat->GetA() << " " << cvol->GetMedium()->GetId() << std::endl;

         std::cout << std::endl; 
         
         l2->SetPoint( i , cpoint[0],cpoint[1],cpoint[2]);

         gGeoManager->FindNextBoundaryAndStep(); // Calculate distance to the next boundary and Evolve system.
         
         i += 1;
  }

   l1->SetLineColor(6);
   l1->SetLineWidth(3);
   l2->SetLineColor(14);
   l2->SetLineWidth(3);
   
   l1->Draw("same");
   l2->Draw("same");


}
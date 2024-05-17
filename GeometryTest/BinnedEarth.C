//Hard Coded
void BinnedEarth()
{


  TCanvas *c1 = new TCanvas("c", "c",0,0,600,600);
  new TGeoManager("simple1", "Simple geometry");
  TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
  TGeoMedium *med = new TGeoMedium("VACUUM",1,mat);
  TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);

  TGeoMaterial *llsvpmat = new TGeoMaterial("EarthMaterial", 1,0.5,1.1*5.51);
  TGeoMedium *LLVP_M = new TGeoMedium("LLSVPComp",1,llsvpmat);



  gGeoManager->SetTopVolume(top);


  std::vector<TGeoMaterial*> MAT; //
  std::vector<TGeoMedium*>  MED; //
  std::vector<TGeoVolume*> LAYER; //

  MAT.push_back(new TGeoMaterial("EarthMaterial", 1,0.5,5.51) );
  MED.push_back(new TGeoMedium("AvgMaterial",1,MAT[0]));

  LAYER.push_back( gGeoManager->MakeSphere("INNERCORE",MED[0], 0,40,0,180,0,360) );//0

  LAYER.push_back( gGeoManager->MakeSphere("OUTERCORE",MED[0], 40,50,0,180,0,360) );//1

  LAYER.push_back( gGeoManager->MakeSphere("LOWERMANTLE",MED[0], 50,60,0,180,0,360) );//2

   //Adding LLSVP
  TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);//Rotation

  TGeoVolume *LLVP = gGeoManager -> MakeSphere("LLVP", LLVP_M , 50, 60 , 0,45,0,360 );

  

  LLVP ->SetLineColor(kGreen);
           
  LAYER[2]->AddNode(LLVP,1,rot1);
  

  LAYER.push_back( gGeoManager->MakeSphere("UPPERMANTLE",MED[0], 60,70,0,180,0,360) );//3

  LAYER.push_back( gGeoManager->MakeSphere("CRUST",MED[0], 70,80,0,180,0,360) );//4
  LAYER.push_back( gGeoManager->MakeSphere("ATMOSPHERE",MED[0], 80,90,0,180,0,360) );//5
  
  //LAYER[0] ->SetLineColor(kRed);
  //LAYER[0] ->SetVisibility(kFALSE);
  LAYER[2] ->SetLineColor(kBlue);
  
  
  LAYER[0] ->SetVisibility(kFALSE);
  //LAYER[1] ->SetVisibility(kFALSE);
  //LAYER[2] ->SetVisibility(kFALSE);
  LAYER[3] ->SetVisibility(kFALSE);
  LAYER[4] ->SetVisibility(kFALSE);
  //LAYER[5] ->SetVisibility(kFALSE);
  

  //LAYER[4] ->SetLineColor(kGreen);



  //Create World
  top->AddNode(LAYER[0],1);
  top->AddNode(LAYER[1],1);
  top->AddNode(LAYER[2],1);
  top->AddNode(LAYER[3],1);
  top->AddNode(LAYER[4],1);
  top->AddNode(LAYER[5],1);



  gGeoManager->CloseGeometry();
  gGeoManager->SetTopVisible(); // the TOP is invisible

  top->Draw();
  TView *view = gPad->GetView();
  view->ShowAxis();
 
 // Tracking----------------------------------------------------------------------------------------------------------------------------------------

  //Calculate Paths inside the Earth

   //Direction of neutrino in spherical coordiates: https://mathworld.wolfram.com/SphericalCoordinates.html
   double zen = 143.0; // zen (90,180]

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

   ROOT::Math::SVector<double, 3> u(xn,yn,zn); // Direction of neutrino
   ROOT::Math::SVector<double, 3> o(0.0,0.0,-90.0); // Origin of the line
   ROOT::Math::SVector<double, 3> cc(0,0,0); // Origin of the Sphere 

   double r = 90;



 
   double a = ROOT::Math::Dot( u, u);
   double b = 2*ROOT::Math::Dot( u, o - cc);
   double c = ROOT::Math::Dot( o-cc, o-cc) - r*r ;
   double d1 = (-b+sqrt(b*b -4*a*c))/(2*a);
   double d2 = (-b-sqrt(b*b -4*a*c))/(2*a);

   std::cout << " d1 " << d1 << " d2 " << d2 << std::endl;

   TPolyLine3D *l1 = new TPolyLine3D();
   TPolyLine3D *l2 = new TPolyLine3D();
   l1->SetPoint( 0 , 0.0, 0.0, -90);
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
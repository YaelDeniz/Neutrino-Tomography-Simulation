//Hard Coded
void BinnedEarth()
{
   /*
   TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   new TGeoManager("sphere", "poza7");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
   gGeoManager->SetTopVolume(top);

   int phibins = 10; //Phi [0 pi]
   int thbins  = 20; //Theta [0 2pi]
   //bin size

   double wth  = (360 - 0)/thbins;
   double wphi = (180 - 0)/phibins;

   double phimin, phimax, thmin, thmax;

   for (int i = 0 ; i < thbins; ++i)
   {
      for (int j = 0; j < phibins; ++j)
      {
         
         phimin = (j)*wphi;
         phimax = (j+1)*wphi;
         thmin  = (i)*wth;
         thmax  = (i+1)*wth;
         TGeoVolume *vol = gGeoManager->MakeSphere("SPHERE",med, 30,40,phimin,phimax,thmin,thmax);
         vol->SetLineWidth(1);
         top->AddNode(vol,1);
         
         std::cout << i << " " << j << std::endl;
      }

   }

  // TGeoVolume *vol = gGeoManager->MakeSphere("SPHERE",med, 30,40,0,180,0,360);
  // vol->SetLineWidth(1);
  // top->AddNode(vol,1);

   //Bin-Layer

   //Make-Layer

   gGeoManager->CloseGeometry();
   //top->Raytrace();
   //gGeoManager->SetNsegments(1);
   top->Draw();
   TView *view = gPad->GetView();
   view->ShowAxis();
   */

  //Define layers 

  TCanvas *c = new TCanvas("c", "c",0,0,600,600);
   
  new TGeoManager("simple1", "Simple geometry");

  TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
  TGeoMedium *med = new TGeoMedium("MED",1,mat);
  TGeoVolume *top = gGeoManager->MakeBox("TOP",med,100,100,100);
  gGeoManager->SetTopVolume(top);

  std::vector<TGeoVolume*> LAYER; //RWH

  LAYER.push_back( gGeoManager->MakeSphere("SPHERE1",med, 30,40,0,180,0,360) );//0
  LAYER.push_back( gGeoManager->MakeSphere("SPHERE2",med, 40,50,0,180,0,360) );//1
  LAYER.push_back( gGeoManager->MakeSphere("SPHERE3",med, 50,60,0,180,0,360) );//2
  LAYER.push_back( gGeoManager->MakeSphere("SPHERE4",med, 60,70,0,180,0,360) );//3
  LAYER.push_back( gGeoManager->MakeSphere("SPHERE5",med, 70,80,0,180,0,360) );//4
  LAYER.push_back( gGeoManager->MakeSphere("SPHERE6",med, 80,90,0,180,0,360) );//5
  
  //LAYER[0] ->SetLineColor(kRed);
  //LAYER[0] ->SetVisibility(kFALSE);
  LAYER[2] ->SetLineColor(kBlue);
  
  
  LAYER[0] ->SetVisibility(kFALSE);
  //LAYER[1] ->SetVisibility(kFALSE);
  //LAYER[2] ->SetVisibility(kFALSE);
  LAYER[3] ->SetVisibility(kFALSE);
  LAYER[4] ->SetVisibility(kFALSE);
  LAYER[5] ->SetVisibility(kFALSE);
  

  //LAYER[4] ->SetLineColor(kGreen);

  //Adding LLSVP
  TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);//Rotation

  TGeoVolume *LLVP1 = gGeoManager -> MakeSphere("LLVP1", med , 53, 55 , 0,45,0,360 );
  TGeoVolume *LLVP2 = gGeoManager -> MakeSphere("LLVP2", med , 60, 68 , 0,25,0,360 );
  

  LLVP1 ->SetLineColor(kGreen);
  LLVP2 ->SetLineColor(kRed);
                
  LAYER[2]->AddNode(LLVP1,1);
  

  //Create World
  top->AddNode(LAYER[0],1);
  top->AddNode(LAYER[1],1);
  top->AddNode(LAYER[2],1);
  top->AddNode(LAYER[3],1);
  top->AddNode(LAYER[4],1);
  top->AddNode(LAYER[5],1);


  gGeoManager->CloseGeometry();

  top->Draw();
  TView *view = gPad->GetView();
  view->ShowAxis();




   


}
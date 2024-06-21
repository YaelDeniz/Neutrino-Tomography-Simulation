/*
In this code:

PremMatrix = [R, RHO, Z/A, LAYER] 
PathMatrix = [L, RHO, Z/A, LAYER] 

R --- Column vector of layer radius
L --- Column vector of Neutrino Baselines


*/



//#include <vector> 
std::vector< std::vector<double> > SortPREMdata( std::string PREM_MODEL = "prem_44layers.txt" )
{
   std::string PREM_PATH = "/home/dehy0499/OscProb/PremTables/"+PREM_MODEL;

   //--- open PREM model 
   std::vector< std::vector<double> > PremMatrix; // Matrix data from Prem model
   std::vector<double> PremData(4);               // [R, rho, z/a, layer]

   // Variables for storing table rows
   double radius, density, zoa, layer;

   std::ifstream PREM_DATA;

   PREM_DATA.open(PREM_PATH);

      // Loop over table rows
   while(PREM_DATA >> radius >> density >> zoa >> layer)
   {

      PremData[0] = radius;
      PremData[1] = density;
      PremData[2] = zoa;
      PremData[3] = layer;

      PremMatrix.push_back(PremData);
    }

   return PremMatrix;

}


void PREM3D( std::vector< std::vector<double> > PremMatrix, Bool_t LLVP = kFALSE, double h_llvp = 1000, double drho_llvp = 1.03, double LLVP_angle = 45.0 )
{

   
   int last = PremMatrix.size()-1; // Index of last Elment

   double Ro = PremMatrix.back()[0]; //Outermost Radius in the Models

   //--Define some transformation
   TGeoTranslation *tr1 = new TGeoTranslation("tr1", 0, 0, 1.0);

   TGeoRotation   *rot1 = new TGeoRotation("rot1", 90.0, 90.0, 0.0);//Rotation
   
   TGeoCombiTrans *combi1 = new TGeoCombiTrans("combi1",0, 0, 1.0, rot1 );
   
   TCanvas *c1 = new TCanvas();
   
   TGLViewer *view = (TGLViewer*)gPad->GetViewer3D();
   
   TGeoManager *geom = new TGeoManager(); // Instance of object 

   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0); // World Material
 
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",9999, matVacuum); // World Medium
    
   TGeoVolume *world = geom -> MakeBox("TOP", Vacuum , 2, 2, 3);

   geom->SetTopVolume(world); //Make top medium the world
     

   //Inner part of Inner Core

   double radmin = 0.0;

   double radmax = PremMatrix[0][0]/Ro;

   double density = PremMatrix[0][1]; 
   
   double zoa = PremMatrix[0][2];
   
   double layer = PremMatrix[0][3];
  
  
    TGeoMaterial *matlay = new TGeoMaterial("matlay",1.0 ,zoa,density);
    //--- define some media
    TGeoMedium *medlay= new TGeoMedium("Layer Material",1.0, matlay);
 


   TGeoVolume *modlayer= geom -> MakeSphere("layer", medlay ,radmin,radmax);

   modlayer->SetVisibility(kFALSE);
   
   world->AddNode(modlayer,99999,tr1);
   
   for (int i = 1; i < PremMatrix.size(); ++i)
   {
      

      radmin = PremMatrix[i-1][0]/Ro;
      
      radmax = PremMatrix[i][0]/Ro;

      double radmid = (radmin + radmax)/2.0; 

      density = PremMatrix[i][1]; 
      
      zoa = PremMatrix[i][2];
      
      layer = PremMatrix[i][3];

      TGeoMaterial *matlay = new TGeoMaterial("mat",1.0 ,zoa,density);
       
      TGeoMedium *medlay = new TGeoMedium("layer Material",layer, matlay);

      TGeoVolume *modlayer = geom -> MakeSphere("layer", medlay , radmin, radmax);

      //std::cout << " radmin " << PremMatrix[i-1][0] << " radmax " << PremMatrix[i][0] << " " << radmid << " R " << Ro <<std::endl;

      modlayer->SetLineColor(kOrange);

      modlayer->SetVisibility(kFALSE);
      
      // Insert LLVPs Sectors

      //int j = 1; //LLVP_index

       if ( LLVP )
       { 
          double local_rho = drho_llvp*density; //Fixed density

          double llvp_rmax = 3480.0 + h_llvp; //Fixed density

          if (radmin*Ro >= 3480.0 && radmax*Ro <= llvp_rmax)
          {

          std::cout << "****************************-*-*-*-*-*" << radmin*Ro << " " << radmax*Ro << " " << (radmax-radmin)*Ro << " " << h_llvp/3 <<  std::endl;
          
             
             TGeoMaterial *matblob= new TGeoMaterial("matblob",1.0 ,zoa,local_rho);
                
             TGeoMedium *medblob = new TGeoMedium("LLVPs Material",7, matblob);
               
             TGeoVolume *LLVP = geom -> MakeSphere("LLVP", medblob , radmin, radmax , 0,LLVP_angle,0,360 );
                
             modlayer->AddNode(LLVP,i+1,rot1);

          }

       }


       world->AddNode(modlayer,i+1,tr1);
   
   }

   geom->CloseGeometry();

   world->SetLineColor(kBlack);

   gGeoManager->SetTopVisible(); // the TOP is invisible

   world->Draw();

}


void Segments()
{
   
   std::vector< std::vector<double> > PathMatrix;
  
   std::vector<double> PathData(4); // Row [L, rho, z/a, layer]
  
   std::vector< std::vector<double> > PremMatrix = SortPREMdata("prem_15layers.txt"); // Sort PREM model into a readable matrix

   double sumr = 0;

   double h = 1000;

   double dh = h/3;

   double Rmax1 = 3480 + dh;
   double Rmax2 = 3480 + 2*dh;
   double Rmax3 = 3480 + 3*dh;

   //double r = 0;

   int k = 0; // Index tracker

   double rmin1,rmax1,rmin2,rmax2,rmin3,rmax3;

   double lb,lt,mb,mt,ub,ut;

   lb = 3480;
   lt = 3480 + dh;

   mt = 3480 + 2*dh;

   ut = 3480 + 3*dh;

   std::vector <double> r;
   //std::vector <double> exces;

//Segment 1-----------------------------------------------------------------------------------------------------------------------
/*
   std::cout << "Segement 1: " << std::endl;

   for (int i = 0; i < PremMatrix.size(); ++i)
   {
      double rc = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (rc >= lb && rc<= lt) //Between lower bottom or lower top
      {
         double rc_top = PremMatrix[i+1][0]; //Upper boundary of current layer

         sumr = sumr + (rc_top - rc); // Current thicknes (sum of layer thickness) 

         r.push_back(dh-sumr); // Reminder height of llvp segment


         //std::cout << "k: " << k << " " <<   rmin << " " << rmax << " " << rmax - rmin << " " << sumr << " " << r[k] << std::endl;


         //If Current thickness is larger than the segment thickness, add the reminder to the current lower boundary.
         if (sumr > dh)
         {
            //sumr = sumr + (h-sumr); 

            std::cout << "top sublayer of segment 1: " << rc << " " << rc + r[k-1] << " " << rc + r[k-1] - 3480 << std::endl;

            mb = rc + r[k-1]; //Create the lower end for the second LLVP segment

            r.clear();

         }

         else
         {
            //std::cout << "sublayer of segment 1: " << k << " " <<   rmin << " " << rmax << " " << rmax - rmin << " " << sumr << " " << r[k] << std::endl;
            std::cout << "sublayer of segment 1: " << rc << " " << rc_top  << std::endl;

         }

         ++k;
         
      }


   }
   std::cout << " " << std::endl;

   std::cout << "Next segment should start at" << mb << std::endl;

   */
// Segment 1 test

   k =0; // restart loop

   double lslabL = 0; //upper slab length

   double srmin = 0;
   double srmax = 0;

   std::cout << " " << std::endl;
   std::cout << " Segement 1: " << lb << " " << lt << std::endl;

   double lbs = lb;

      for (int i = 0; i < PremMatrix.size(); ++i)
      {

      double layerb = PremMatrix[i][0]; //Lower boundary of current layer

      
       
      if (layerb > lb )
      {
         std::cout << "*Current layer " << layerb << std::endl;


         lslabL = lslabL + (layerb - lbs );

         

         if (lslabL > dh)
         { 
            double excess = dh - lslabL;
            srmin= lbs;
            srmax= layerb + excess;
            std::cout << "top of lower segment: " << srmin << " " << srmax << " " << lslabL + excess << " " << srmax-srmin << std::endl;
            mb = srmax;
            r.clear();
            break;
         }
    
         else
         {
            srmin = lbs;
            srmax = layerb;
            r.push_back(dh-lslabL);
            std::cout << srmin << " " << srmax << " " << lslabL << " " << srmax-srmin << std::endl;
         }
         lbs = layerb;
        
        
         
         ++k;
  
      }


   }


//Segment 2------------------------------------------------------------------------------------------------------------------------
   k =0; // restart loop

   double mslabL = 0; //mid slab length

   std::cout << "Segement 2: " << mb << " " << mt << std::endl;

   
   double mbs = mb;

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double layerm = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (layerm > mb )
      {


         mslabL = mslabL + (layerm - mbs );

         //r.push_back(dh-mslabL);

         if (mslabL > dh)
         {
            double excess2 = dh-mslabL;
            srmin= mbs;
            srmax= layerm + excess2;
            std::cout << "top of middle segment: " << srmin << " " << srmax << " " << mslabL + excess2 << " " << srmax-srmin << std::endl;
            ub = srmax;
            r.clear();
            break;
         }

         else
         {
      

            srmin = mbs;
            srmax = layerm;
            r.push_back(dh-lslabL);
            std::cout << srmin << " " << srmax << " " << mslabL << " " << srmax-srmin << std::endl;

         }
         mbs = layerm;
        
        
         
         ++k;




         
      }


   }

   //Segment 3------------------------------------------------------------------------------------------------------------------------
   k =0; // restart loop

   double uslabL = 0; //mid slab length

   std::cout << "Segement 3: " << ub << " " << ut << std::endl;



   double ubs = ub;

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double layeru = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (layeru > ub )
      {


         uslabL = uslabL + (layeru - ubs );

         //r.push_back(dh-mslabL);

         if (uslabL > dh)
         {
            double excess3 = dh-uslabL;
            srmin= ubs;
            srmax= layeru + excess3;
            std::cout << "top of upper segment: " << srmin << " " << srmax << " " << uslabL + excess3 << " " << srmax-srmin << std::endl;
            r.clear();
            break;
         }

         else
         {
      

            srmin = ubs;
            srmax = layeru;
            r.push_back(dh-uslabL);
            std::cout << srmin << " " << srmax << " " << uslabL << " " << srmax-srmin << std::endl;

         }
         ubs = layeru;
        
        
         
         ++k;




         
      }


   }
/*
//Segement 3--------------------------------------------------------------------------------------------------------------------

   k =0; // restart loop

   double uslabL = 0; //upper slab length

    srmin = 0;
    srmax = 0;

   std::cout << " " << std::endl;
   std::cout << " Segement 3: " << ub << " " << ut << std::endl;

   double ubs = ub;

      for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double lub = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (lub >= ub )
      {


         uslabL = uslabL + (lub - ubs );

         

         if (uslabL > dh)
         {
            srmin= ubs;
            srmax= ubs + r[k-1];
            std::cout << "top of upper segment: " << srmin << " " << srmax << " "<< " " << uslabL << " " << srmax - srmin << " " <<  r[k] << std::endl;
            r.clear();
            break;
         }

         else
         {
            srmin = ubs;
            srmax = lub;
            r.push_back(dh-uslabL);
         std::cout << srmin << " " << srmax << " " << uslabL << " " << srmax-srmin << " " <<  r[k] << std::endl;
         }
         ubs = lub;
        
        
         
         ++k;




         
      }


   }

*/  
 }



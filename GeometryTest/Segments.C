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
  
   std::vector< std::vector<double> > PremMatrix = SortPREMdata("prem_44layers.txt"); // Sort PREM model into a readable matrix

   double sumr = 0;

   double h = 1000;

   double dh = h/3; //Slab Thickness

   double PileThickness = 1000;

   double thicknessOfSection = PileThickness/3.0; // Salabas have same thickness

   double Rmax1 = 3480 + dh;
   double Rmax2 = 3480 + 2*dh;
   double Rmax3 = 3480 + 3*dh;

   //double r = 0;


   double rmin1,rmax1,rmin2,rmax2,rmin3,rmax3;

   double lb,lt,mb,mt,ub,ut;

   double lowerBottom, lowerTop, midBottom, midTop, upperBottom, upperTop;

   lowerBottom = 3480;

   lowerTop = 3480 + dh;

   midTop = 3480 + 2*dh;

   upperTop = 3480 + 3*dh;

// Segment 1 test


   double slabLength =0;

   double InnerR = 0;
   double OuterR = 0;

   std::cout << " " << std::endl;
   std::cout << " Segement 1 radius: " << lowerBottom << " - " << lowerTop << std::endl;

   double lowerSegment = lowerBottom; //Lower Radius of a segment in Pile section

      for (int i = 0; i < PremMatrix.size(); ++i)
      {

      double RadiusInBottom = PremMatrix[i][0]; //Lower boundary of current layer

      
       
      if (RadiusInBottom > lowerBottom )
      {
         std::cout << " Radius of layer above lower bottom  " << RadiusInBottom << std::endl;


         slabLength = slabLength + (RadiusInBottom - lowerSegment ); //Calculate segement lenght 

         if (slabLength > thicknessOfSection)
         { 
            double excess = thicknessOfSection - slabLength;
            InnerR= lowerSegment;
            OuterR= RadiusInBottom + excess;
            std::cout << "top of lower segment: " << InnerR << " " << OuterR << std::endl;
            midBottom = OuterR;
            break;
         }
    
         else
         {
            InnerR = lowerSegment;
            OuterR = RadiusInBottom;
            //r.push_back(dh-lslabL);
       std::cout << "Radius of lower section's segements: "<< InnerR << " " << OuterR << std::endl;
         }

         lowerSegment = RadiusInBottom;
        
      }


   }


//Segment 2------------------------------------------------------------------------------------------------------------------------
   slabLength =0;
   InnerR = 0;
   OuterR = 0;
   
   std::cout << " " << std::endl;
   std::cout << " Segement 2 radius: " << midBottom << " - " << midTop << std::endl;

   
   double midSegment = midBottom; //Lower Radius of a segment in Pile section

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInMiddle = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (RadiusInMiddle > midBottom )
      {
         std::cout << " Radius of layer above middle bottom  " << RadiusInMiddle << std::endl;


         slabLength = slabLength + (RadiusInMiddle - midSegment ); //Calculate segment lenght 


         //r.push_back(dh-mslabL);

         if (slabLength > thicknessOfSection)
         {
            double midExcess = thicknessOfSection-slabLength;
            InnerR= midSegment;
            OuterR= RadiusInMiddle + midExcess;
            std::cout << "top of middle segment: " << InnerR << " " << OuterR << std::endl;
            upperBottom = OuterR;
            break;
         }

         else
         {
      

            InnerR = midSegment;
            OuterR = RadiusInMiddle;
       std::cout << "Radius of middle section's segements: "<< InnerR << " " << OuterR << std::endl;

         }

         midSegment = RadiusInMiddle;

      }

   }

   //Segment 3------------------------------------------------------------------------------------------------------------------------
   slabLength =0;
   InnerR = 0;
   OuterR = 0;
   
   std::cout << " " << std::endl;
   std::cout << " Segement 3 radius: " << upperBottom << " - " << upperTop << std::endl;



   double upperSegment = upperBottom;

   for (int i = 0; i < PremMatrix.size(); ++i)
   {

      double RadiusInTop = PremMatrix[i][0]; //Lower boundary of current layer
       
      if (RadiusInTop > upperBottom )
      {

         std::cout << " Radius of layer above uppper bottom  " << RadiusInTop << std::endl;

         slabLength = slabLength + (RadiusInTop - upperSegment );

         //r.push_back(dh-mslabL);

         if (slabLength > thicknessOfSection)
         {
            double topExcess = thicknessOfSection-slabLength;
            InnerR= upperSegment;
            OuterR= RadiusInTop + topExcess;
            std::cout << "top of upper segment: " << InnerR << " " << OuterR << std::endl;
            break;
         }

         else
         {
            InnerR = upperSegment;
            OuterR = RadiusInTop;
            std::cout << "Radius of top section's segements: "<< InnerR << " " << OuterR << std::endl;
         }
         upperSegment = RadiusInTop;
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



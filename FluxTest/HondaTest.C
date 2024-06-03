void FluxTest()
{
    int xbins = 1000;
    int ybins = 1000;
    TCanvas *c1 = new TCanvas();
    TH2F *h = new TH2F("h","h",xbins,0.,10.,ybins,0.,10.);
    TF2 *xyg = new TF2("xyg","xygaus",0.0,10,0.0,10);
    xyg->SetParameters(1,5,2,5,2); //amplitude, meanx,sigmax,meany,sigmay 
    h->FillRandom("xyg",xbins*ybins); 
    h->Draw();

    

    //TF2 *f2 = new TF2("f2","sin(x)*sin(y)/(x*y)",0,1,0 ,1);

   // hist->FillRandom("f2",400000);

    /*
    TRandom *rand =  new TRandom(10);

    for (int i = 0; i < 1e6 ; i++)
    {

        double x = rand->Gaus();
        double y = rand->Gaus();

        hist->Fill(x,y);

    }

    TH2F *h = new TH2F("h","h",100,0.,10.,100,0.,10.);
    TF2 *xyg = new TF2("xyg","xygaus",0,10,0,10);
    xyg->SetParameters(1,5,2,5,2); //amplitude, meanx,sigmax,meany,sigmay 
    h->FillRandom("xyg"); 
    h->Draw();





    */ 
   
   TRandom *rand =  new TRandom(10);

    double xrad = rand->Uniform(10);
    double yrad = rand->Uniform(10);

 

    std::cout << xrad << " " << yrad <<  std::endl;

    std::cout << "True val: " <<  xyg->Eval(xrad,yrad) << " Interpolated: " << h->Interpolate(xrad,yrad) <<  std::endl;

	
    gStyle->SetPalette(kRainBow);
    
    h->SetContour(1000);
    
    h->Draw("COLZ");

    TCanvas *c2 = new TCanvas();
    xyg->Draw();
    
    c1->Print("./Hist2DFunc.png");
    c2->Print("./funcTEST.png");
   



}

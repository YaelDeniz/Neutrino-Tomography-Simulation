void FluxTest()
{

    TCanvas *c1 = new TCanvas();
    TH2D *hist = new TH2D("hist","Histogram",100,0,1,100,0,1);
    hist->SetStats(0);

    TF2 *f2 = new TF2("f2","sin(x)*sin(y)/(x*y)",0,1,0 ,1);

    hist->FillRandom("f2",400000);

    /*
    TRandom *rand =  new TRandom(10);

    for (int i = 0; i < 1e6 ; i++)
    {

        double x = rand->Gaus();
        double y = rand->Gaus();

        hist->Fill(x,y);

    }
    */
   
   TRandom *rand =  new TRandom(10);

    //double xrad = rand->Uniform();
    //double yrad = rand->Uniform();

    double xrad = 0.5;
    double yrad = 0.5;

    std::cout << xrad << " " << yrad <<  std::endl;

    std::cout << "True val: " <<  f2->Eval(xrad,yrad) << " Interpolated: " << hist->Interpolate(xrad,yrad) <<  std::endl;

	
    gStyle->SetPalette(kRainBow);
    
    hist->SetContour(1000);
    
    hist->Draw("COLZ");

    TCanvas *c2 = new TCanvas();
    f2->Draw();
    
    c1->Print("./Hist2DFunc.png");
    c2->Print("./funcTEST.png");
   



}

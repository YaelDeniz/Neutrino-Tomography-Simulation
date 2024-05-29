void FluxTest()
{

    TCanvas *c1 = new TCanvas();
    TH2D *hist = new TH2D("hist","Histogram",100,-1,1,100,-1,1);
    hist->SetStats(0);

    TRandom *rand =  new TRandom(10);

    for (int i = 0; i < 1e6 ; i++)
    {

        double x = rand->Gaus();
        double y = rand->Gaus();

        hist->Fill(x,y);

    }

    hist->Draw();

    hist->Draw("COLZ");
    c1->Print("./Hist2D.png");
   



}

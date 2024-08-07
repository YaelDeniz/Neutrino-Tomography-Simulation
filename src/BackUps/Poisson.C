/*
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

using namespace std;


int main()
{


    std::random_device rd;
    std::mt19937 gen(rd());
 
    // If an event occurs 4 times a minute on average, how
    // often is it that it occurs n times in one minute?
    std::poisson_distribution<> d(2.5);
 
    for (int n = 0; n != 10; ++n)
        std::cout << d(gen) << std::endl; 
        
    
    return 0;

}
*/



#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"


void Poisson()
{   

    double mu = 4.0;

    std::random_device rd; //produce non-deterministic random numbers.

    std::mt19937 gen(rd()); //high quality, but not cryptographically secure, unsigned integer random numbers of type UIntType on the interval
 
    // If an event occurs 4 times a minute on average, how
    // often is it that it occurs n times in one minute?
    std::poisson_distribution<> d(mu);

    TH1F *hist =  new TH1F("hist","histogram",50,0,50);

    //TF1 *fit = new TF1("fit","gaus", 0 , 50);

    TF1 *f1 = new TF1("f1","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, 50); // "xmin" = 0, "xmax" = 10
    
    f1->SetParameters(1, 1, 1); // you MUST set non-zero initial values for parameters

    for (int i = 0; i < 10000; ++i)
    {
        std::cout << d(gen) << std::endl;
        hist->Fill(d(gen));

    }

    TCanvas *c1 = new TCanvas();
    //Double_t factor = 1.;
    //hist->Scale(factor/hist->Integral());
    hist->Draw();
    hist->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1";

}
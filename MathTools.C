/* Implemantation file*/

#include "MathTools.h"

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <valarray>
#include <string.h>

// Calculates the factorial of n
double factorial(const int& n)
{
	double f = 1;
	int i=0;
	for(i = 1; i <= n; i++)
	{
		f *= i;
	}
	return f;
}

// Generate a random number from Poisson distribution
double  Sample_Pois(double lambda, double u)
{ 

	//std::cout << "Integration" << std::endl;

	const int N = 100; // Max Number of Occurances

	//Create CDF for Poisson distribution 

	double prob=0,tmpProb=0,fact=0;

	/* List of variables

	"prob" -> Stores the result of evaluating cumulative probability of the poisson distribution at certain value : prob[k] = Fx(K < k)
	"tmpProb" -> An object that only is used to construct the CDF.
	"fact" -> Result of Factorial k!
	"u" -> Is a random number between 0 to 1 extracted from Unif(0,1)*/

	int sample=0; //Random Sample from Poisson Distribution.
	
	double p[N+1]; //Stores CDF

	// probability of no events
	p[0] = exp(-lambda);
	prob = p[0]; 

	//CMF for Poisson distribution
	for (int k=1; k<N+1; k++)
	{

		fact = factorial (k);
		tmpProb = pow(lambda, k)*exp(-lambda) / fact;
		prob += tmpProb;
		p[k] = prob;
		
	}

	//Generate Sample from image using CDF( F[x] ): F[ x < m-1]<= u  < F[x < m]   
	for (int m=100; m>=0; m--)
	{
			 
		if ( u < p[m] )
		{	

			sample = m; 
				
		}

	}

	return sample; 
}

// Integrals using Newton Cotes
double ncquad(std::vector<double> xdat, std::vector<double> ydat){

    double h = xdat[1]-xdat[0];
    int m = (xdat.size() - 1)/2 ; 
    int n = 2*m; //Number of intervals
    double f[n + 1],x[n + 1]; // Value of f at x_i

    for(int k=0; k < n + 1 ; ++k)
    {

        x[k]=xdat[k];
        f[k]=ydat[k];
    
    }

	// Numerical Integration Composite Simpson's Rule

    double I_f = 0;
    I_f += f[0];

    for (int j = 1; j <= m ; ++j)
    {

        I_f += 4*f[2*j-1]; //odd index

    }
    for (int j = 1; j <= m-1 ; ++j)
    {

           I_f += 2*f[2*j];

    }

    I_f += f[n+1];;
    I_f = (h/3)*I_f;
    return I_f;
}

double simpson(std::vector<double> x, std::vector<double> f) //Simpsons Rule 1/3

{
  

  double I = 0;

  double h = (x[2]-x[0])/2;

  I = (h/3)*( f[0] + 4*f[1] + f[2] ); // Simpsons Rule 1/3

  return I;


}

// Produce even space points in LogSpace (inputs are exponents).

std::vector<double> logspace(double a, double b, int k) {

  const auto exp_scale = (b - a) / (k - 1);
  std::vector<double> logspace;
  logspace.reserve(k);
  for (int i = 0; i < k; i++) {
    logspace.push_back(i * exp_scale);
  }
  std::for_each(logspace.begin(), logspace.end(),
                [](double &x) { x = pow(10, x); });
  return logspace;
}
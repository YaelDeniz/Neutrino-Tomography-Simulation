
#ifndef MATHTOOLS_H
#define MATHTOOLS_H

/* 

This files contains declarations and provides Interface to other 
parts of a program 

*/


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <valarray>


//Factorial of a function
double factorial(const int& n);

//Generates random Samples from Poisson Distribution
double  Sample_Pois(double lambda, double u);

//Numericla Integration
double ncquad(std::vector<double> xdat, std::vector<double> ydat);

//Logarithm Spaced points (Even spaced points in LogSpace)
std::vector<double> logspace(double a, double b, int k);

#endif

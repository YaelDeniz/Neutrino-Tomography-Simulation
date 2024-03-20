#include "PhyTools.h"

// Neutrino Xsec approximations
double XSec(double E, int ptype)
{           

    double xsc = 0;

    double c = 1E-4; // Conversion factor cm2 -> m2

    switch(ptype)
    {
        //case "neutrino":
        case 1:
                xsc = (0.75E-38*E)*c; // Results is in m2 
                break;

        //case "antineutrino":
        case -1:
         //       std::cout << "NUBAR" << std::endl;
                xsc = (0.35E-38*E)*c; // Results is in m2 
                break;
    }

    return xsc;
}
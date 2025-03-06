void AdjustableBins(){

	double ErecoMin = 2;
	double ErecoMax = 10;
	
        std::vector<double> ErecoEdges; //cos(theta) binning
        ErecoEdges.push_back(ErecoMin);
        double Erecoi = ErecoMin;


	double Ae = 0.05;
	double Be = 0.1;
	
	std::cout << "Finder" << std::endl;
	
	

	while (ErecoEdges.back() < ErecoMax) {
	
	
	
		//Find center
		TF1 *frecoC = new TF1("fmin","-1*x + [0] + 0.5*( [1]*x+[2]*sqrt(x) )", 0, 100);
		frecoC->SetParameter(0,Erecoi);
		frecoC->SetParameter(1,Ae);
		frecoC->SetParameter(2,Be);
		ROOT::Math::RootFinder finder_center;
	    	finder_center.Solve(*frecoC ,0,100);
	    	double ErecoC = finder_center.Root();
		
		//Find next edge
		double Erecoii = ErecoC + 0.5*(Ae*ErecoC + Be*sqrt(ErecoC));
		
		ErecoEdges.push_back(Erecoii);
		
		Erecoi = Erecoii;
		
		std::cout << "Value finder" << Erecoi << std::endl;
		
		}
		
		double ErecobinEdges [ErecoEdges.size()];
		
		for (int i=0 ; i < ErecoEdges.size(); i++ )
		{
		
		ErecobinEdges[i] = ErecoEdges[i];
		
		std::cout << "Bin edge" << ErecobinEdges[i] << std::endl;
		
		}
		
		//Create Histogram 
		int nrecoBins = sizeof(ErecobinEdges)/sizeof(ErecobinEdges[0]) - 1; // Number of bins

	    	// Create a histogram with variable bin widths
	    	TH1D *h = new TH1D("h", "Histogram with Variable Bin Width", nrecoBins, ErecobinEdges);
	    	
	    	for (int j = 1; j <= h->GetNbinsX(); j++) 
	    	{
		std::cout << "Bin " << j << ": [" << h->GetBinLowEdge(j) 
		      << ", " << h->GetBinLowEdge(j+1) << "," << h->GetBinLowEdge(j+1) - h->GetBinLowEdge(j) << "," << h->GetBinWidth(j) << "]" << std::endl;
	       
		}
		
		
	       std::cout << "End finder" << std::endl;

}


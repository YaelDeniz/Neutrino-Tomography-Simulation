void string2double()
{

        // Input string containing doubles
    std::string input = "3.14 2.71 1.61 0.577 1.41";
    
    // String stream for parsing the input string
    std::istringstream iss(input);
    
    // Vector to store the parsed doubles
    std::vector<double> doubles;
    
    // Temporary variable to store each parsed double
    double temp;
    
    // Parse the input string and store the doubles in the vector
    while (iss >> temp) {
        doubles.push_back(temp);
    }
    
    // Example of assigning parsed doubles to individual variables
    if (doubles.size() >= 5) {
        double d1 = doubles[0];
        double d2 = doubles[1];
        double d3 = doubles[2];
        double d4 = doubles[3];
        double d5 = doubles[4];

        // Print the variables to verify the result
        std::cout << "d1: " << d1 << "\n";
        std::cout << "d2: " << d2 << "\n";
        std::cout << "d3: " << d3 << "\n";
        std::cout << "d4: " << d4 << "\n";
        std::cout << "d5: " << d5 << "\n";
    } 
    else 
    {
        std::cerr << "Not enough doubles in the input string." << std::endl;
    }

}
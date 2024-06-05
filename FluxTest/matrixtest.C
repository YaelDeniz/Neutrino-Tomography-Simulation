std::vector<std::vector<int>> reshape(const std::vector<int>& vec, int a, int b) {
    // Check if the total number of elements matches
    if (vec.size() != a * b) {
        throw std::invalid_argument("The size of the vector does not match the dimensions of the matrix.");
    }

    // Initialize the reshaped matrix
    std::vector<std::vector<int>> matrix(a, std::vector<int>(b));

    // Fill the matrix with elements from the vector
    for (int i = 0; i < a; ++i) {
        for (int j = 0; j < b; ++j) {
            matrix[i][j] = vec[i * b + j]; //Row Major Order
        }
    }

    return matrix;
}

void matrixtest() {
    // Example input vector
    std::vector<int> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    int a = 3;
    int b = 4;

    // Reshape the vector into a matrix
    try {
        std::vector<std::vector<int>> matrix = reshape(vec, a, b);

        // Print the reshaped matrix
        for (const auto& row : matrix) {
            for (const auto& element : row) {
                std::cout << element << " ";
            }
            std::cout << std::endl;
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
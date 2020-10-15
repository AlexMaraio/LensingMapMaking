#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

extern "C" {

    void alm_to_cl(const char* file_path, double* clm)
    {
        // Open the provided file
        std::ifstream input_file(file_path);

        // Raise an error if the file is not open
        if(!input_file.is_open()) { throw std::runtime_error("Cannot open the specified file! Is the path correct?"); }

        // Temp variables for data to be read into
        unsigned int ell, m;
        double Re, Im;       

        // Go through each line individually
        while (input_file >> ell >> m >> Re >> Im)
        {
            // Compute the modulus-square of the alm
            auto alm_sq = Re * Re + Im * Im;

            // If m is not zero, we need to multiply by two for the complex conjugate
            if (m != 0) {alm_sq *= 2;}

            // Add this alm to the Cl list with normalisation
            clm[ell] += alm_sq / (2 * ell + 1);
        }
        
        // Close the file once we're done with it
        input_file.close();
    }

}
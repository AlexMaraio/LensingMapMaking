#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <complex>
#include <chrono>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <omp.h>

#include "ComputeAlms.h"


// Functions visible to Python need to be C compatible
extern "C" {

void alm_to_cl(const char* file_path, double* clm)
{
  // Open the provided file
  std::ifstream input_file(file_path);

  // Raise an error if the file is not open
  if(!input_file.is_open())
  { throw std::runtime_error("Cannot open the specified file! Is the path correct?"); }

  // Temp variables for data to be read into
  unsigned int ell, m;
  double Re, Im;

  // Go through each line individually
  while(input_file >> ell >> m >> Re >> Im)
  {
    // Compute the modulus-square of the alm
    auto alm_sq = Re * Re + Im * Im;

    // If m is not zero, we need to multiply by two for the complex conjugate
    if(m != 0) { alm_sq *= 2; }

    // Add this alm to the Cl list with normalisation
    clm[ell] += alm_sq / (2 * ell + 1);
  }

  // Close the file once we're done with it
  input_file.close();
}

void map_to_alm(const char* file_path)
{
  // Get the number of threads available using OpenMP
  const int num_threads = omp_get_max_threads();
  std::cout << "Number of threads: " << num_threads << "\n";

  // Open the input file provided
  std::ifstream input_file(file_path);

  // Since this file had a header which we don't need, read the first line to a dummy string
  std::string dummy_line;
  getline(input_file, dummy_line);

  // The input file has six columns, and so we need six temporary variables
  double theta, phi;
  double gamma1_z1, gamma2_z1;
  double gamma1_z2, gamma2_z2;

  // Create six vectors which the data will be stored in
  std::vector<double> thetas, phis;
  std::vector<double> gamma1_z1s, gamma2_z1s;
  std::vector<double> gamma1_z2s, gamma2_z2s;

  constexpr unsigned int ell = 76;
  constexpr unsigned int m = 75;

  // Constants associated with the specific Flask run. TODO: make these dynamic.
  constexpr auto n_side = 2048;
  constexpr auto n_pix = 12 * n_side * n_side;

  std::cout << "Reading in data\n";
  // Read in the data and store them in the given vectors
  while(input_file >> theta >> phi >> gamma1_z1 >> gamma2_z1 >> gamma1_z2 >> gamma2_z2)
  {
    thetas.emplace_back(theta);
    phis.emplace_back(phi);
    gamma1_z1s.emplace_back(gamma1_z1);
    gamma2_z1s.emplace_back(gamma2_z1);
    gamma1_z2s.emplace_back(gamma1_z2);
    gamma2_z2s.emplace_back(gamma2_z2);
  }

  std::cout << "Data read in\n";

  //? Single-threaded approach
  auto alm = std::complex<double>(0);

  auto start_time = std::chrono::system_clock::now();

  // Go through each datapoint and compute the value of Y_lm(theta, phi) * f(theta, phi) and add to alm
  for(std::size_t i = 0; i < thetas.size(); ++i)
  {
    alm += std::conj(boost::math::spherical_harmonic(ell, m, thetas[i], phis[i])) * gamma1_z1s[i];
  }

  // Normalise the alm value through the area of each pixel
  alm *= (4 * M_PI) / n_pix;

  std::cout << "alm value using boost: " << alm << "\n";

  double elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::system_clock::now() - start_time).count();
  std::cout << "First way took: " << elapsed_time << "\n";


  //? Single-threaded approach but now using my own class for speed comparison

  // Initiate Alm class
  auto AlmClass = ComputeAlms(105, n_side, thetas, phis, gamma1_z1s);

  // Time evaluation of a single alm using my class
  start_time = std::chrono::system_clock::now();
  std::cout << "alm from class: " << AlmClass.calc_alm(ell, m) << "\n";
  elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::system_clock::now() - start_time).count();
  std::cout << "Second time took: " << elapsed_time << "\n";

  //? Compute a list of alms:
  start_time = std::chrono::system_clock::now();
  auto alms = AlmClass.calc_alms();
  elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::system_clock::now() - start_time).count();

  for(const auto& alm_i: alms) { std::cout << alm_i << "\n"; }
  std::cout << "Alms computation time took: " << elapsed_time / 1000 << " seconds, with an average time of "
            << elapsed_time / alms.size() << " milliseconds per alm \n";

  AlmClass.output_data("/home/amaraio/Documents/PhD/Codes/LensingMapMaking/Data/AlmOutput.dat");
}

}


int main()
{
  std::cout << "Minimum value of custom type: " << boost::math::tools::min_value<my_mpfr_float>() << "\n";
  std::cout << "Maximum value of custom type: " << boost::math::tools::max_value<my_mpfr_float>() << "\n";

  map_to_alm("/home/amaraio/Documents/PhD/Codes/LensingMapMaking/Data/Non-linear/Output-shear-map.dat");

  return 0;
}

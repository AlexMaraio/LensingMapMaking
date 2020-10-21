//
// Created by amaraio on 19/10/2020.
//

#ifndef COMPUTEALMS_H
#define COMPUTEALMS_H

#include <stdexcept>
#include <cmath>
#include <vector>
#include <complex>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/multiprecision/mpfr.hpp>


// Definition of our custom Boost multiprecision type which will be used in calculations
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<15>> my_mpfr_float;

// * In order to use boost::multiprecision, we need to overload the sqrt operator.
namespace std
{
    using boost::multiprecision::sqrt;
}

// Function to obtain the index of an (ell, m) pair in a 1D flattened list
// Also used to get the number of (ell, m) configurations up to ell_max
inline auto idx(const int ell, const int m)
{
  return ell * (ell + 1) / 2 + m;
}


// Function that returns the square-root prefix for a spherical harmonic Ylm
template<class number=long double>
inline auto harmonic_prefix(const unsigned ell, const unsigned m)
{
  // The prefix for the spherical harmonics have two factorial terms, which we compute here using the Gamma function.
  auto ratio = boost::math::tgamma<number>(ell - m + 1) / boost::math::tgamma<number>(ell + m + 1);

  ratio *= (2 * ell + 1) / (4 * M_PI);

  ratio = std::sqrt(ratio);

  return ratio;
}

// Function that returns the complex phase factor of a spherical harmonic as a complex number x + iy.
// Note: We have a negative imaginary part as we are interested in the conjugate of the spherical harmonic.
inline auto complex_exp(const int m, const double phi)
{
  return std::complex<double>(std::cos(m * phi), -std::sin(m * phi));
}

// Function that wraps the Boost associated Legendre polynomial function.
// Note: return type is of the multiprecision variable
template<class number=long double>
inline auto get_legendre_poly(const int ell, const int m, const double theta)
{
  return boost::math::legendre_p<number>(ell, m, std::cos(theta));
}


// Class that stores all of the information, and does the computation of an alm calculation
class ComputeAlms
{
public:

    // Constructor
    ComputeAlms(const int ell_max, const int n_side, const std::vector<double>& thetas, const std::vector<double>& phis,
                const std::vector<double>& in_data) :
            ell_max(ell_max),
            N_side(n_side),
            N_pix(12 * n_side * n_side),
            thetas(thetas),
            phis(phis),
            input_data(in_data),
            alms(std::vector<std::complex<double>>(idx(ell_max, ell_max) + 1))
    {
      // Check that the input arrays are all the same size
      if((thetas.size() != phis.size()) || (thetas.size() != in_data.size()))
        throw std::runtime_error("The input arrays must all have the same size");
    }

public:
    // Compute a single Ylm function
    inline auto compute_Ylm(int ell, int m, double theta, double phi);

    // Compute an alm value
    template<class number=long double>
    auto calc_alm(int ell, int m);

    // Calculate a set of alm values up to ell_max
    template<class number=long double>
    std::vector<std::complex<double>> calc_alms();

    // Save the alm data to a file given by file_path
    void output_data(std::string file_path) const;


private:
    // Input data
    const int ell_max;
    const int N_side;
    const int N_pix;

    // Vectors storing the input data
    const std::vector<double>& thetas;
    const std::vector<double>& phis;
    const std::vector<double>& input_data;

    // Vector which we will store the result in
    std::vector<std::complex<double>> alms;
};

// Generic function to return the value of a spherical harmonic Y_lm evaluated at (theta, phi)
auto ComputeAlms::compute_Ylm(const int ell, const int m, const double theta, const double phi)
{
  const auto Ylm = harmonic_prefix(ell, m) * get_legendre_poly(ell, m, theta);
  return std::complex<double>(static_cast<double>(Ylm * std::cos(m * phi)),
                              static_cast<double>(-Ylm * std::sin(m * phi)));
  // return harmonic_prefix(ell, m) * get_legendre_poly(ell, m, theta); // * complex_exp(m, phi);
}

// Calculate the alm value for an (ell, m) configuration
template<class number>
auto ComputeAlms::calc_alm(const int ell, const int m)
{
  // Create our alm variable which we will add to
  number alm_re = 0;
  number alm_im = 0;

  // Variables for storing the previous theta value, and hence the value of the Legendre polynomial at that theta
  double prev_theta = 0;
  number prev_legendre = 0;

  // Go through each data point in the provided data arrays
  for(std::size_t i = 0; i < thetas.size(); ++i)
  {
    // If our new theta is not the same as the previous one, store this along with the Legendre polynomial at this theta
    if(prev_theta != thetas[i])
    {
      prev_theta = thetas[i];
      prev_legendre = get_legendre_poly<number>(ell, m, prev_theta);
    }

    // Add this value to the alm
    alm_re += prev_legendre * std::cos(m * phis[i]) * input_data[i];
    alm_im -= prev_legendre * std::sin(m * phis[i]) * input_data[i];
  }

  // Normalise through the area of each pixel
  alm_re *= (4 * M_PI) / N_pix;
  alm_im *= (4 * M_PI) / N_pix;

  // Multiply our alm by the alm's prefix
  alm_re *= harmonic_prefix(ell, m);
  alm_im *= harmonic_prefix(ell, m);

  // Return the alm value casting it to a double as we no longer need the increased precision.
  return std::complex<double>(static_cast<double>(alm_re), static_cast<double>(alm_im));
}

// Function that computes all alm values up to ell_max using OpenMP parallelization
template<class number>
std::vector<std::complex<double>> ComputeAlms::calc_alms()
{
  for(int ell = 80; ell <= ell_max; ++ell) //! Changed this to start at 80
  {
    #pragma omp parallel for // NOLINT(openmp-use-default-none)
    for(int m = 0; m <= ell; ++m)
    {
      alms[idx(ell, m)] = calc_alm<number>(ell, m);
    }
  }

  return alms;
}

void ComputeAlms::output_data(const std::string file_path) const
{
  // Create the output stream object
  std::ofstream out_file(file_path, std::ios_base::out | std::ios_base::trunc);

  // Write the header
  out_file << "ell" << "\t" << "m" << "\t" << "Re" << "\t" << "Im" << "\n";

  // Go through each (ell, m) combination and print it's value
  for(int ell = 80; ell <= ell_max; ++ell) //! Again, changed this
  {
    for(int m = 0; m <= ell; ++m)
    {
      const auto alm = alms[idx(ell, m)];
      out_file << ell << "\t" << m << "\t" << alm.real() << "\t" << alm.imag() << "\n";
    }
  }

  // Close the file once done
  out_file.close();
}

#endif //COMPUTEALMS_H

# Quadratic Maximum Likelihood estimators

## Estimators with spatially-varying noise levels

When we have added noise to our shear maps we have always assumed that the variance of the noise is constant
with position. However, this assumption may not necessarily be correct, as the noise levels can change over the sky,
and so we want to include this in our analysis. To model this at a beginners level, we simply induce a multiplicative
bias correction into the noise standard deviation when we generate the noise: we multiply it by 0.8 for the northern
hemisphere, and by 1.2 for the southern hemisphere. One realisation of this noise map looks as follows

![Spatially varying noise map](figures/Eclipse/VaryingNoise/NoiseMap_Q.png)

Here, we can see that the amplitude of the shape noise is generally larger in the southern hemisphere over the northern
one.

With our new noise maps in place, we can look at how this affects the recovered power spectra between the Pseudo-Cl
and QML methods. We can average over 256 maps to find the means and variances as

![Varying noise avg cls](figures/Eclipse/VaryingNoise/Cl_avg_ratio_11.png)

![Varying noise var cls](figures/Eclipse/VaryingNoise/Cl_var_ratio_11.png)

Here, we are also comparing the averages and variances between the spatially-varying noise levels and a constant noise
level. This shows that both estimators are capable of recovering the input power spectrum, on average, despite having
varying noise levels. We also see that for our chosen bias scales, there is little difference in the variances
between the regular noise and varying noise cases. It was noted, however, that when the bias terms became much larger,
then the QML estimator significantly out-performed Pseudo-Cl, and so that could be something to look into.  
Looking at the averages and variances, it is unlikely that the spatially-varying noise levels (at least with our bias
values), will have any meaningful impact on recovered cosmological parameters.

## Parameter constraints in the As-ns plane

With our power spectrum estimation techniques working well, we can now apply them to an ensemble of maps to generate
parameter constraint contours. Here, a Gaussian likelihood has been used where we have neglected any correlations
between different _l_ values. Results are presented for maps with an N_side of 64, so an l_max of 128, for two redshift
bins.

![As histogram](figures/Eclipse/Paramerter/As_Histogram.png)

![ns histogram](figures/Eclipse/Paramerter/ns_Histogram.png)

![As-ns KDE](figures/Eclipse/Paramerter/Asns_KDE.png)

Here, we see very little difference between QMl and PCl method in either the 1D histograms or the 2D combined
distribution of the mean parameter values. It would be worth, however, applying our estimates to higher-dimensional
problems where we have many more parameters to estimate and marginalise over at the same time.

## Testing the LSST's Core Cosmological Library

So far, I have been using CAMB to compute the theory spectra for the shear maps at specific redshifts and as a function
of cosmological parameters in the MCMC analyses. However, this is quite a slow process taking multiple seconds to
evaluate the Cl values per run, and so for a general MCMC analysis this needs to be sped up. One way to do this is to
use a cosmological emulator, which are many orders of magnitude faster than calling CAMB every single time. One such
emulator is the LSST's Core Cosmological Library (CCL) code. This can evaluate the shear power spectra much faster than
CAMB, ~9 seconds for CAMB compared to ~75 ms for CCL, and should provide just as accurate values. We can compare the
two codes by taking the ratio of the Cl values for what should be the same cosmology and redshift bin combinations:

![Ratio of CCL to CAMB](figures/CCL_CAMB_ratio.png)

Here, we see that there are quite large differences at the lowest _l_ values, and also quite a difference with the 
high-_l_ values too.

## How to downgrade maps effectively

We know that the QML method only works for low resolution maps (N_side of 64 is about its maximum), and so some form
of combined estimator of QML and PCl could be used. Hence, we need a way to downgrade our shear maps from the normal
high resolution to ones that can be used in our QML analysis. So far, we have been using HealPix's `ud_grade` function
which simply averages over the smaller pixels to form the larger pixel's value. This is not ideal as it suppresses power
on small scales, which we would like to keep in our maps:

![Diff N_side Cls](figures/Cl_NsideAvg_ratio.png)

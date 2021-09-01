# Fisher Matrices

So far, when trying to get parameter constraints I have been using an ensemble of MCMC simulations and averaging their results.
This is quite a slow process as each MCMC simulation can take a while to run, with this increasing dramatically when
many parameters are being sampled over, and the need for many simulations to be run to get accurate ensemble statistics.
Instead, we can use the parameter's Fisher matrix to find approximate parameter constraints much faster as we don't need
to run costly MCMC simulations and can simply evaluate the Fisher matrix. Here, we are evaluating the Fisher matrix
as written in Equation 128 of [1910.09273](https://arxiv.org/pdf/1910.09273.pdf), which evaluates the derivatives of the
Cl coefficients for each parameter for each redshift bin combination and each ell mode independently.

To ensure that we have coded up the Fisher matrix correctly, we can compare the parameter constraints obtained from
the Fisher matrix to those from an MCMC run, which gives

![MCMC vs Fisher for 1024](figures/Fisher/FisherMCMCComparison_1024.png)

Note that we have used the covariance matrix for the MCMC simulation, and so does not necessarily capture the
non-Gaussian nature of the MCMC analysis. Here, we see excellent agreement between the two methods, and so I am confident
in my Fisher matrix implementation.

## Using Fisher matrices to predict 3x2pt data

The above implementation of the Fisher matrix was for three redshift bins using cosmic shear only. This gave us six
unique spectra, and so while the Cl covariance matrix was a 6x6 matrix, it was possible to compute and code-up by hand.
However, as we want to extend our analysis to galaxy counts spectra and the galaxy-shear cross-spectra, this would give
us 21 unique spectra. This would give us a 21x21 Cl covariance matrix, which is certainly unfeasible to compute by-hand.
Thus, an automated method to compute it was needed and implemented efficiently.  
Using this, we can compute the individual Fisher matrices for the cosmic shear spectra and then galaxy number counts
spectra only, and then the combined Fisher matrix. This gives us

![Individual vs combined 3x2pt comparison](figures/Fisher/Fisher3x2ptComparison.png)

Here, we see that cosmic shear only provides good constraints on Omega_c and _h_, whereas number counts only provide
good constraints on sigma_8, and so combining both probes (and then adding their cross-spectra) provides the tightest
constraints. It is also interesting to see that the degeneracy direction differs between the probes, which is
another useful thing when combining probes.

## Applying this to dark energy

Now that we have a general formalism to compute the Fisher matrix, we can apply this to any combination of cosmological
parameters that `CCL` accepts. Here, we can use it to look at how the dark energy constraints as a function of the
number of ell-modes in our Cl data vector. When we consider three values of ell_max, of 128, 512, and 2048, we find

![Dark energy constraints using 3x2pt](figures/Fisher/Fisher_3x2pt_w0wa.png)

We can then turn these contours into a single value called the Figure-of-Merit, which encodes how tight the w0-wa
contour is, and gives

![Dark energy figure-of-merit](figures/Fisher/DarkEnergyFoM.png)

This shows that the FoM increases with some power-law as a function of the maximum ell-mode available, and so for any
analysis it is important to include as many ell-modes as possible.

## Pseudo-Cl vs QML for parameter constraints

Here, we are taking the Cl Fisher matrix for the QML method and the Gaussian covariance matrices for the Pseudo-Cl
method as the Cl covariance matrices for our two different methods and running them through my parameter Fisher pipeline.
Here, we are using both the QML method where we have used the (very costly) analytic Fisher matrix and the approximation
of the Fisher matrix using the covariance matrix of the y<sub>_l_</sub> values.

![Fisher matrix comparison for PCl vs QML](figures/Fisher/Fisher_3x2pt_PClQML_32.png)

## Properly implementing the Cl Fisher matrix

In my previous QML calculations, I used a naive assumption that the Cl values for a spectra XY was composed of the 
y_l values for that spectra deconvolved with the specific Fisher matrix elements belonging to XY. This is not the case,
as the Cl values for any spectra are a result of the deconvolution of the _entire_ array of y_l values with the **full**
Fisher matrix. My approximation works in the noise-dominated case, as seen before, as then the off-diagonals in the
Fisher matrix are suppressed. With the full pipeline now in place, the results for the Cl variances are as follows:

### TT

![Var of the TT](figures/ShearTEB/TT_VarRatio.png)

### EE

![Var of the EE](figures/ShearTEB/EE_VarRatio.png)

### TE

![Var of the TE](figures/ShearTEB/TE_VarRatio.png)

Ignoring the right-hand edge of all these plots, where we are above _l_ > 2 * N<sub>side</sub>, we see that the ratio 
of the variances for all three non-zero spectra are consistently above one, which indicates that QML might not tend to
the Pseudo-Cl estimator at high _l_ values, as previously thought. Hence, it is important to push our QML estimators to
larger and larger resolutions to see if this trend continues, or the ratio decays to one eventually.

## Propagating the correct variances to constraints on dark energy

Now that my QML pipeline is working correctly for one redshift bin, we can compute the Fisher matrix for the dark
energy parameters up to an _l_<sub>max</sub> of 128, giving us

![Dark energy Fisher, with correct vars](figures/Fisher/Fisher_3x2pt_TEB_Masked.png)

Here we see moderately decreased contours for the QML method over Pseudo-Cls, which we would expect from looking at
the variances alone. However, it is to be noted that the Fisher matrix calculation uses the variances and covariances
of all combinations of Cls, and thus even though the variances of the Cl values may be smaller for QML over PCl,
the off-diagonal terms in the Cl covariance matrix still need to be studied.

## Comparisons of the Figure-of-Merit for different spectra and N<sub>side</sub>

Now that I have a pipeline that can produce estimates of the dark energy Figure-of-Merit (FoM) for different spectra
(TT only, EE only, TT EE and TE combined), we can investigate how the FoM changes with resolution for these different 
probes:

![Dark energy figure-of-merit for different spectra](figures/Fisher/DarkEnergyFoM_Nside.png)

Here we see some interesting behaviour: on a log-log plot it appears that the FoM increases linearly with the
resolution of the map, with QML being slightly more constraining than PCl - though the difference shrinks with 
resolution. We also see that initially, the TT spectra were more constraining than the EE modes, however they both
converg to the same accuracy when increasing the resolution.


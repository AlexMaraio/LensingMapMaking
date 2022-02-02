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

## Fisher forecasts using numerical Cl-Fisher matrices

Now that we have an implementation of our numerical Fisher matrices, we can extend our analysis of how the figure of
merit changes between QML and PCl to larger values of Nside.

### The number of maps doesn't seem to matter

When computing the Fisher matrix via finite-differences, it is essential to average over enough maps to get an accurate
prediction for the Cl Fisher matrix. If we look at using only five maps, we can take the ratio of the numerically-estimated
and the analytically-evaluated Cl-Fisher matrices, to get

![Ratio of Fisher matrices using 5 maps](figures/FiniteDiffFisher/FisherRatio_5maps.png)

Here, we can see that these ratios are quite large, with over a ±60% error, which is not very precise at all! However,
if we increase the number of maps to five hundred (and thus increase the computational time by a factor of one hundred),
we find the accuracy of the obtained Fisher matrix is vastly improved:

![Ratio of Fisher matrices using 500 maps](figures/FiniteDiffFisher/FisherRatio_500maps.png)

Here, we see an order of magnitude decrease in the ratio, which shows that it's much more accurate.

Knowing that the Cl-Fisher matrix's accuracy vastly improves when averaging over more maps, we can ask if this difference
propagates into the parameter's Fisher matrix. Surprisingly, we find that it makes an almost negligible difference here
when comparing 5 and 500 maps:

![Comparison of param-Fisher for 5 vs 500 maps](figures/FiniteDiffFisher/FisherComparison_nummaps.png)

Here, we see that all parameter contours are practically on top of each other, and thus we can use only a handful of
maps to average over.

### Extending Nside

Previously, we were limited by RAM constraints to maps with Nside 64 or less. However, we can extend to at least an
Nside of 128 using our new method, giving:

![Parameter Fisher matrices](figures/FiniteDiffFisher/ParamFisher_QML_vs_PCl.png)

![Figure of merit vs Nside](figures/FiniteDiffFisher/FoM_vs_Nside.png)


## Looking at the diagonals of the Fisher matrix

Above, we have seen that as we increase the number of maps that we average over, the accuracy of numerical Cl-Fisher
matrix improves. Here, we specifically want to see how the diagonal values (which carry most of the signal in the
Fisher matrix) responds to an increase in the number of maps used:

![Diagonal of Cl-Fisher](figures/FiniteDiffFisher/FisherDiag_64.png)

![Ratio of diagonal of Cl-Fisher](figures/FiniteDiffFisher/FisherDiagRatio_64.png)

Here, we can see that the noise gradually decreases with _l_, and indicates that at high-_l_ (where the signal is 
largest), we are more accurate.

## Figure of merit for dark energy

Here, we want to see how the figure of merit (FoM) changes as a function of N_side between QML and PCl using different
apodization scales:

![FoM vs N_side](figures/FiniteDiffFisher/FoM_vs_Nside_w0wa.png)

![FoM vs N_side ratio](figures/FiniteDiffFisher/FoM_vs_Nside_w0wa_ratio.png)

Here, we see that when we apodize the mask we get a roughly constant ratio between PCl and QML.

## Using correct cosmic shear signal

As previous parts of my code used the 'TT' spectra as galaxy number counts for use in a 3x2pt analysis, this was used
when going to a single scalar field map instead of the convergence signal. These two spectra are quite different, as
shown in the plot below:

![Correct and wrong TT signals](figures/FiniteDiffFisher/TT_EE_Cls.png)

Here, we can see that the previously used signal was much larger than the correct signal, which meant that the 
signal-to-noise ratio was much larger than it should have been. Thus, when computing the conjugate-gradient of the 
covariance matrix, the signal-to-noise ratio term was far more dominant than it should have been - and thus taking a 
lot longer to converge.  
With the correct signal in place, the signal-to-noise ratio is much lower and so the covariance matrix becomes more
diagonal, thus easier to invert and so takes less time to converge. This allows us to push the resolution much further
than it was previously able to do; here we go to an Nside of 512 (taking over 15 hours to do so, though) with 1024 being
potentially possible.

![Param Fisher using Nside of 512](figures/FiniteDiffFisher/ParamFisher_N512_combined.png)

At this resolution, we see a dramatic effect of apodizing the mask, even only on a scale of a single degree 
**when using the additional star mask**.

## Extending numerical Fisher estimate to EB shear

Now that I have an accurate way to estimate the Cl-Fisher matrix for our spin-0
TT field, this can be easily extended to compute the Fisher matrix for our spin-2
EE, EB, and BB fields. 

### Comparing numerical values to analytic result

First, we want to check that the extension of my C++ code to deal with the
spin-2 fields is correct, through checking that the numerical values roughly
agree with the analytic result. 

![Ratio of numerical to analytic Fisher matrix](figures/FiniteDiffFisher/Spin2/NumericalFisherFull_AnalyticRatio_N64.png)

Here, we inject power into modes along the x-axis, and then measure the response
in the modes on the y-axis, averaged over five realisations. Note that we can't
inject power into just the EB modes when we have no B-modes, and so their 
numerical values are very close to zero.  
This shows that our EE- and BB-mode estimation are pretty good, as their ratios
are very close to one, with what looks like random noise populating their 
matrices.

### Extending to Nside of 256

Now that we have verified that our code is working, we can extend it to an
Nside of 256, which took about ten hours to compute:

![Full numerical Fisher matrix](figures/FiniteDiffFisher/Spin2/NumericalFisher_N256.png)

Since the results for the EB-modes are very unreliable, as they have a magnitude
that is many orders of magnitude smaller than the other two modes, we can
restrict our Fisher matrix to just EE and BB modes:

![Cut-down numerical Fisher matrix](figures/FiniteDiffFisher/Spin2/NumericalFisherCut_N256.png)

Here, we have also explicitly symmetrised the Fisher matrix by adding its 
transpose and dividing by two, which makes it numerically invertible. 

### Comparing with Pseudo-Cl coupling matrices

We can now compare our inverse Fisher matrix to the Pseudo-Cl coupling matrix 
as computed from NaMaster, which for the full EE & BB matrix gives us 

![Ratio of PCl coupling mat to inverse numerical Fisher matrix](figures/FiniteDiffFisher/Spin2/PCl_invFisher_ratio_N256_E0_B0.png)

This shows that the covariances are generally much larger for the Pseudo-Cl's 
vs our QML ones, especially at low l values and for the EB & BB modes in
particular. To get an idea of how the variances differ between our two methods,
we can simply take the diagonal of our individual covariance ratios, which gives

![Diag-ratio of PCl coupling mat to inverse numerical Fisher matrix](figures/FiniteDiffFisher/Spin2/PCl_invFisher_ratiodiag_N256_E0_B0.png)

Here, we see the very large increase in errors for the EB spectra, and that 
the auto-spectra of EE and BB tend to a constant value, which is due to the 
apodisation of the mask.

#### Purifying B-modes

If we now enable B-mode purification, then we find the Pseudo-Cl coupling matrix
is significantly modified, which gives our ratio & diagonal values as

![Ratio of PCl coupling mat to inverse numerical Fisher matrix w/ B-mode purify](figures/FiniteDiffFisher/Spin2/PCl_invFisher_ratio_N256_E0_B1.png)

![Diag of ratio of PCl coupling mat to inverse numerical Fisher matrix w/ B-mode purify](figures/FiniteDiffFisher/Spin2/PCl_invFisher_ratiodiag_N256_E0_B1.png)

Here, we see much larger errors induced in both the EE and BB spectra when
purifying the B-modes, especially at low-l, and then tend to the same value
(which, again, is due to the requirement of apodisation to use B-mode 
purification). For context, the variance ratio of just the TT spectra is

![Diag of ratio of PCl coupling mat to inverse numerical Fisher matrix - T only](figures/FiniteDiffFisher/Spin2/PCl_invFisher_ratiodiag_N256_T.png)

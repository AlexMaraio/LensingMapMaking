# Effects of dark energy on the lensing power spectrum

One of the main goals of measuring cosmic shear is to make accurate and precise measurements of dark energy, and
hopefully provide constraints on the type of dark energy that we observe.
In order to do this, we need to understand how dark energy affects the theoretical lensing power spectrum.

As weak lensing is sensitive to the evolution of the matter distribution, we want to see how this varies as we
look into different models of dark energy. Here, we compare the standard ΛCDM model, where dark energy is a
cosmological constant with equation-of-state parameter _w_=-1, with the *w*CDM model, where the dark energy
equation-of-state varies with scale factor *a* as
<center>
<img src="https://render.githubusercontent.com/render/math?math=w(a) = w_0 %2B w_a(1-a)," height=20> 
</center>

as suggested in [astro-ph/0208512](https://arxiv.org/pdf/astro-ph/0208512.pdf).

Plotted below is the ratio of the total matter and dark energy densities as a function of scale factor, normalised
to the critical density, for both models where *w*0 = -0.9 and *w*a = 0.1.

![Omega matter and Lambda](figures/Omega_mLambda.png)

Here, we see that the evolution of the dark energy density is higher in the wCDM model over ΛCDM, which corresponds
to a slightly smaller matter density between *a*~0.1 and *a*~1. Since the matter density is smaller, we expect
to find slightly suppressed growth of structure in this model which would lead to a smaller lensing power spectrum.  
The lensing power spectra for these two models is plotted below, and follows what we expect given the above plot.

![Lensing power spec w0 wa](figures/LensingPowerSpec_w0wa.png)

Here, we find that the wCDM model predicts a slightly smaller amplitude of the lensing power spectra, which makes
sense given the smaller matter density evolution in this model. We note that the two models predict the same shape
of the power spectrum, just slightly different amplitudes.

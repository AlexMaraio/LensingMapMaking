# Investigating the effects of masks

A mask is simply where we neglect observations taken over a region of the sky. This is often done because the observations
in this region are heavily contaminated with noise, and so the data is not coming from the underlying signal that we are
trying to detect. An example of this is when we mask over the Milky Way in most observations due to the increased stellar
density, which dramatically increases the noise in the signal.

## Masking CMB data

First, I will look at masks in relation to CMB data as there are many existing demonstrations of masks applied to WMAP and
*Planck* data online. Once the basics are understood, I can apply this to weak lensing data.

Here, we plot raw temperature anisotropy data obtained by WAMP seven-year data in the Q-band (41 GHz), with a normalised colour-map scale.

![Galaxy counts power spectrum](figures/investigating_masks/WMAP_no_mask.png)

We see that the data is heavily contaminated by the Milky Way running through the centre of the image, and so if we want
to do any science with this map, we need to mask over this area. This will leave behind the actual temperature anisotropies
that are from the CMB, and not from any other sources. To do so, WMAP employed the following mask

![Galaxy counts power spectrum](figures/investigating_masks/WMAP_mask.png)

In this figure, a value of one corresponds to allowing the data through the mask and a value of zero corresponds to
the data being masked out. Here, we see the large section that is masked out in the centre is the Milky Way. However,
we also see many smaller regions have also been masked out throughout the entire map, which corresponds to bright stars and
other objects that produce a large amount of noise - and so need to be masked out accordingly.

With this mask in place, we can now plot the raw temperature anisotropy data of the residual signal, which looks like

![Galaxy counts power spectrum](figures/investigating_masks/WMAP_with_mask.png)

Hence, we recover a much cleaner signal with values that match the theoretical prediction of the CMB signal alone.

## Masks for Euclid

Above, we were looking at masks for WMAP, which of course was a CMB experiment and so much less sensitive to bright sources
or nearby galaxies than optical galaxy surveys. We can now look at what a prototype Euclid mask may look like, which is

![Galaxy counts power spectrum](figures/investigating_masks/Euclid_mask.png)

Here, we can see that there is a large band around the centre of the map, which corresponds to masking out the Milky Way,
as this plot is in galactic coordinates. We also see that there is a slightly smaller band that appears to oscillate in our
mask. This is the ecliptic plane of the Solar System, and looks this way because we are plotting the mask in galactic
coordinates. If we were plotting this mask in ecliptic coordinates, then we would find that the Solar System runs
through the centre of our mask, by definition.

From a mask, we can define the fraction of the surveyed sky *f* as

<center>
<img src="https://render.githubusercontent.com/render/math?math=f_\textrm{sky} \equiv \frac{\textrm{Unmasked&#160;area}}{\textrm{Total&#160;area}} ," height=40> 
</center>

For the Euclid trial mask above, we find the fraction of sky observed is around 37%.

### Creating a semi-realistic Euclid mask

As discussed above, the two most important objects when making a mask for an optical galaxy survey are the Milky Way and the
Solar System. Hence, we can create a semi-realistic mask by just masking them out in their respective coordinate systems.  
This is shown below:

![Milky way and solar system mask](figures/investigating_masks/MilkyWay_SolarSystem_mask.png)

Here, we have taken the Milky Way to be a cylindrical band between θ = π/3 and θ = 2π/3 in galactic coordinates, and the
Solar System to be another cylindrical band but one slightly thinner going from θ = 2π/5 to θ = 3π/5 in ecliptic coordinates.

> Do note that these are just arbitrary angles for both objects that I have plucked out of thin air based on the rough size
> of the Euclid map above!

It is important to ask why we need to mask out both these regions. For the Milky Way, we have many bright star sources
which affects the quality of the signal, but importantly for both the Milky Way and the Solar System, there are large
quantities of cosmic dust in both of these regions which attenuates the observed signal. This means that we will observe
less galaxies in these regions, and so the data analysis becomes much harder to perform as we don't know what the "true"
signal should be in these situations.

Now that we have our two individual regions that we want to mask out, we need to combine these to form a single mask.
To do so, we need to set up a HealPy `Rotator` object to rotate our Solar System mask from ecliptic to galactic coordinates.
Now that we have both masks in a single coordinate system, we can combine them to form a new mask which only allows data
through if the region is unmasked in both maps. As these maps are simply NumPy arrays, we can use `np.logical_and` to
combine them in the correct way. This gives us our combined map, which in galactic coordinates, looks like

![Combined dummy mask](figures/investigating_masks/Combined_dummy_mask.png)

Here, we see our dummy mask looks very much like the prototype Euclid mask (which makes sense, and is good!) except for a
few extra regions which our very simple approximation of cylindrical regions does not include.  
We also note that our fraction of sky observed using the dummy mask is around 32%, which is broadly comparable to the
Euclid mask.

## Obtaining power spectra from masked maps

Now that we have a realistic mask for Euclid-like observations, we can now apply this to maps and see what their effects
are on the power spectrum coefficients recovered from the map.

First, let us apply the Euclid mask to some noisy convergence data provided from Flask.

![Combined dummy mask](figures/investigating_masks/Converg_with_mask.png)

Now that we have our masked map, we can compute what the *Cl* values are using the masked map

![Combined dummy mask](figures/investigating_masks/Converg_power_with_mask.png)

Here, we have plotted the *Cl*s recovered from the **full** map in blue, and those recovered from the masked map in green.
Nothing that the masked *Cl* values look like the original *Cl* values, just at a lower amplitude, we have applied a very
naive "correction" to the masked *Cl*s where we simply divide by the fraction of sky observed, which results in the orange
line. We see that this basic correction seems to work remarkably well across all *l* values, subject to a bit of noise,
which is very interesting to see.  

Now that we have some recovered *Cl* values from our mask (with our basic correction in place), we can look at how these
recovered values compare to the original, unmasked values. Here, we tale the ratio of the masked values with the unmasked
values, to give the plot below.

![Combined dummy mask](figures/investigating_masks/Ratio_of_converg_power.png)

Here we see that the recovered *Cl* values are quite noisy, with respect to the original values, especially at low *l*.
This somewhat dissipates at higher *l*, but still quite large differences.

## Power spectra for an ensemble of masked maps

Above, we have looked at the power spectrum for a single map that had a mask applied. To get more detailed information
about how applying a mask changes the distribution of the *Cl* values obtained from a masked map, we need to repeat this
process many times to build up an ensemble of masked maps, and then can compare summary statistics between the masked
and unmasked distributions.

### Mean, variance, skew, and kurtosis

First, we look at how the first four moments of our *Cl* distribution change when we apply a mask. This is for data
that has been collected for one thousand runs in each dataset.  
Note that throughout this section we renormalise the recovered *Cl* values from the masked map by dividing by the fraction
of sky let through.

![Combined dummy mask](figures/investigating_masks/ensemble_avg/Average.png)

Here we are plotting the deviation of the average of the recovered *Cl* values with respect to the input *Cl* values.
Interestingly, we see that for *l* below a hundred, the masked *Cl*'s on average predict a slightly lower power spectrum,
with this difference decreasing at larger *l*.

![Combined dummy mask](figures/investigating_masks/ensemble_avg/Variance.png)

Here we are plotting the variance of the recovered set of *Cl* values divided by the squared-average. We see that the
masked values predict a systematically larger variance than the unmasked values, which follow the cosmic variance
prediction almost exactly.

![Combined dummy mask](figures/investigating_masks/ensemble_avg/Skew.png)
![Combined dummy mask](figures/investigating_masks/ensemble_avg/Kurtosis.png)

Here, we see that the skew and kurtosis statistics are generally much noisier than the mean or variance, but we see
a general trend that the masked values predict a slightly larger value for both.

We can try and get a better visualisation of the increased variance by producing KDE plots at each *l* for the two
sets of values.

![Combined dummy mask](figures/investigating_masks/ensemble_avg/Ridge_plot.png)

### Looking at the covariance values

Now that we have our set of masked *Cl* values, we can look at if using a mask introduces any correlation between the
*Cl* values at different *l*. Previously, we have found that for unmasked maps, there is very little correlation between
the *Cl* values, with a maximum correlation coefficient of |r| ~ 0.06. However, we can repeat the same analysis but now
for our masked values, which gives:

![Combined dummy mask](figures/investigating_masks/ensemble_avg/Covariance.png)

Here, we see that there is a very interesting peak that lies on the diagonal between *l* and *l* - 2. Here, we see a strong
positive correlation between these values that was not previously observed for the unmasked values.

**Note:** I have updated this plot with a new colour-bar which better reflects the nature of this data. Here, we still
see the extremely strong correlation between *l* and *l* ± 2 modes, and again very little correlation outside of these
modes.

## Looking at *Cl* bleed between *l*'s for masked maps

Above, we can see that for our masked maps, there seems to be a slight positive correlation between different *l* modes
for our masked maps. We can investigate this further by creating maps where we start our with just a single *Cl* value
non-zero, applying a mask and then recovering the *Cl* values from this.  
To do this, we start off with two maps where all *Cl* values are zero, expect for *l* of 50 and 250. The maps for these
to *Cl* values are shown below.

![Combined dummy mask](figures/investigating_masks/Cl_bleed/Dummy_maps.png)

Now that we have these two maps, we can apply the Euclid mask and recover what the power spectrum is for the masked
maps, which is shown below.

![Combined dummy mask](figures/investigating_masks/Cl_bleed/Recov_power_spec.png)

## Power spectra for an ensemble for different masks

Here, we want to look at the statistical features in the power spectrum by applying different sized masks. To do so, we
need to average over vary many different realisations of individually masked maps to be able to view the affects of the
mask, and isolate statistical noise to the best that we can.

### Constructing nine masks

First, we needed to construct a variety of masks that have different sky fractions. To do so, we generalised the previous
method of constructing galactic and ecliptic cuts for the Milky Way and Solar System, respectively, but dynamically
changing the thickness of the two cuts. Implementing this over a range of different thicknesses gives us the following
nine masks

![Nine masks](figures/investigating_masks/many_masks/9masks.png)

We can see that our process has generated masks that range from blocking nearly all of the pixels through to letting
most of the pixels through. By applying these different masks to many different convergence maps, we can try to find
any statistical trends in the *Cl*'s on the sky-fraction.

### Editing Flask

In order to implement the masking and recovery of the *Cl* values for our many different masks in the most efficient way
possible, the source code of Flask was edited to include the desired functionality.  
The result of this was that Flask would recover eleven sets of *Cl* values for each realisation of the convergence map
(the nine masks above, plus the Euclid mask, and of course the unmasked values), and save each to the disk. These were
then read in by the Python code to be added to existing DataFrames.  
With N_side of 1024 and l_max of 1000, both for performance reasons, the code took an average of 5.8 seconds per sample
to do this.

### Results

We now present the results of averaging over our ensemble of maps, where we have generated 6,000 maps.

#### Deviation from input power spectra

Here, we want to look at how the average power spectra recovered from the masked maps deviates from the input power
spectrum. Note that we have used a rolling average using three *l* modes at each point, which helped smooth out some
of the rapid oscillations present in the recovered data, leaving behind the underlying pattern.

![Nine masks](figures/investigating_masks/many_masks/RelativeDifference.png)

Here, we see some very interesting behaviour: for low-*l* modes (less than around 200) we see that there is reduced power
for the masked maps compared to the original power spectrum. We see that the more heavily the map is masked (indicated by
a lower f<sub>sky</sub>), the more the power is suppressed on these large angular scales.
Then above this point, there seems to be slightly more power compared to the original power spectrum,
again with this increase growing as we decrease f<sub>sky</sub>.

#### Variance

We now want to look at how the variance of the recovered *Cl*'s change as we increase the masking of the maps. This
is shown in the following figure

![Nine masks](figures/investigating_masks/many_masks/Variance.png)

Here, we have plotted the expectation from the Γ-distribution as a dashed cyan line. Here we see that as we decrease
f<sub>sky</sub>, the variance seems to increase with a constant factor across all *l*. We can quantify this increase, by
taking the ratio of the recovered lines with the expected values, to give

![Nine masks](figures/investigating_masks/many_masks/VarianceRatio.png)

Here, we see that for all lines except for the most heavily masked maps, this ratio is basically constant across all *l*,
and is some function of f<sub>sky</sub>. A few functions were tried, e.g. 1/f<sub>sky</sub> and 1/sqrt(f<sub>sky</sub>),
but nothing tired *yet* seems to fit these values.

#### Skew

We now want to look at the third moment of the *Cl*'s, skew. Again, we have taken a rolling average of three *l* modes
to damp down some of the oscillations present in the original data.

![Nine masks](figures/investigating_masks/many_masks/Skew.png)

Here, we see that as we decrease f<sub>sky</sub>, the skewness increases, showing that even at very large *l* the
skewness deviates significantly from the Gaussian case and from the expected values from the Γ-distribution (which
is shown in the dashed cyan line).

#### Kurtosis

We can now look at the third moment of the *Cl*'s, kurtosis. Again, we have taken a rolling average of three *l* modes
to clean up the data somewhat.

![Nine masks](figures/investigating_masks/many_masks/Kurtosis.png)

Again, we see that as we decrease f<sub>sky</sub>, the kurtosis of the data increases. However, once *l* increases past
around 100, all lines except the most heavily masked, reduce to basically zero kurtosis (except noise).

### Signal-to-noise of masked maps

A very simple signal-to-noise calculation was performed for the nine masks above, averaging over our ensemble of
six thousand realisations. Here, we have used the theoretical covariance matrix from the Γ-distribution, which assumes
that all *l*-modes are independent of each other. This is true for our unmasked map, but we note that modes start to
mix with each other as we introduce a mask. The results of this are

![Signal to noise](figures/investigating_masks/Signal-to-noise.png)

Here, we see that as we increase the fraction of sky surveyed, the S/N ratio increases, which is what we expected.

## Recovering cosmological parameters from masked maps

Now that we have established methods to mask convergence maps and naively recover the unmasked power spectrum, we can
apply these methods to recover the cosmological parameters that were used to generate the original map. This is a
formulation of a **very basic** likelihood code, which seeks to evaluate the posterior as a function of a singular
cosmological parameter to find the maximum-likelihood value.  

First, though, we need to define what the likelihood is in our case. Here, we have one set of "observed" *Cl*'s, denoted
with a hat, and want to find the most likely parameters of our model that go into generating our theory *Cl*'s. If we
assume that the *Cl*'s are normally distributed (we have shown this not to be true, but we will apply this approximation
anyway), then we find that for a single *Cl* value, the logarithm of the posterior probability is proportional to

<center>
<img src="https://render.githubusercontent.com/render/math?math=-\ln f(C_\ell \,|\, \hat{C}_\ell) \propto (2 \ell %2B 1) \left[\ln(C_\ell) %2B \frac{\hat{C}_\ell}{C_\ell} \right]." height=45> 
</center>

Since the total likelihood is proportional to the multiplication of *f* over all *l*, the log-likelihood will simply
be the sum of the above quantity over all *l*:

<center>
<img src="https://render.githubusercontent.com/render/math?math=-\ln \mathscr{L} \propto \sum_{\ell} (2 \ell %2B 1) \left[\ln(C_\ell) %2B \frac{\hat{C}_\ell}{C_\ell} \right]." height=45> 
</center>

Here, we have assumed a uniform prior on the theory *Cl*'s (and thus cosmological values), but this can be incorporated
into the above equation if so desired.

We can now implement the above equation into our code, and extract results from it.

### Recovering A_s

The first parameter that was attempted to be recovered from a map was the scalar amplitude A_s. Here, we computed the
likelihood for the unmasked map, the Euclid-like mask, and our first mask with f<sub>sky</sub> = 0.15%. This allows us
to see if applying a mask changes the recovered cosmological parameters, even for our very simple case.

The first test was on a pre-existing convergence map which was generated with A_s = 2.1E-9, so we want to see how
close we can get to this value using our maximum-likelihood technique.

![Nine masks](figures/recovering_params/A_s1.png)

Here, we can see that the maximum-likelihood values predicted by all three maps are very close to the original value
(which is shown by the dashed cyan line in the zoomed insert). We can see that the affect of masking gradually
decreases the recovered A_s parameter, with the results being

* Unmasked map: A_s = 2.102E-9
* Euclid mask: A_s = 2.099E-9
* Crazy mask: A_s = 2.094E-9

Hence, even our crazy mask predicts a value of A_s that is very close to the true value, which is encouraging!

We can now repeat the same analysis, but for a true A_s of 2.25E-9, just to see how robust the predictions are.

![Nine masks](figures/recovering_params/A_s2.png)

Here, we now see that the unmasked and Euclid masked maps predict A_s values that are very close to the true value,
with our crazy mask predicting somewhat of a lower value of A_s, but still in the right ballpark:

* Unmasked map: A_s = 2.252E-9
* Euclid mask: A_s = 2.251E-9
* Crazy mask: A_s = 2.215E-9

### Recovering the neutrino masses

Above, we were just trying to recover A_s, which simply scales the entire power spectrum and so is quite a simple value
to predict from a map. Here, we wish to extend the above analysis, but now try to recover the sum of the neutrino masses.
This causes scale-dependant suppression in the power spectrum, and so may be harder to predict than a simple scaling.

Again, we look at the unmasked map, the Euclid-like mask, and our crazy mask.

![Nine masks](figures/recovering_params/m_nu1.png)

Here, we see some interesting results: all likelihoods dip around m_nu = 0.02eV, which causes the maximum-likelihood
value for the unmasked map to be much smaller than what it should be. The two masked maps also feature this dip,
but with is slightly weaker, and so the maximum-likelihood values here are much closer to the input value of 0.06eV.
The results for this are:

* Unmasked map:  sum m_nu = 0.00214 eV
* Euclid mask: sum m_nu = 0.0587 eV
* Crazy mask:  sum m_nu = 0.0620 eV

Hence, even with a much harder value to predict from a map, our maximum-likelihood technique still manages to recover
sensible, and quite close, values for the two masked maps.

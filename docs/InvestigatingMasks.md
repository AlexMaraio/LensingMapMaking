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

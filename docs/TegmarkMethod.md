# Tegmark method for recovering Cl's

Here, we have computed a simple implementation of the Tegmartk method to recover 
the Cl values for a masked map without loosing any information, as described in 
[`astro-ph/9611174`](https://arxiv.org/pdf/astro-ph/9611174.pdf). 

The results of this are

![Unmasked map](figures/Tegmark/TegmarkUnmasked.png)

![Unmasked map](figures/Tegmark/TegmarkMasked.png)

Here, we can see that the masked power spectrum is exactly the same as the original unmasked
power spectrum, which shows that the method is indeed lossless.

If we take the average of the power spectrum for 25 maps, then we find

![Averaged map](figures/Tegmark/MultipleRuns.png)



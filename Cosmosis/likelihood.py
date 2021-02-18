from cosmosis.datablock import names, option_section

import numpy as np
from scipy import constants as sciconst
from scipy import interpolate as sciinterp
from camb.sources import GaussianSourceWindow
import healpy as hp
import camb
from camb.dark_energy import DarkEnergyPPF

# We have a collection of commonly used pre-defined block section names.
cosmo = names.cosmological_parameters


def setup(options):
    # The lmax value that we want to generate the Cl values up to
    lmax = 2000
    ells = np.arange(2, lmax + 1)

    # The N_side parameter of the generated maps & masks
    n_side = 2048

    # Convert our N_side parameter to the number of pixels in the map
    n_pix = 12 * n_side ** 2

    use_noise = options[option_section, "use_noise"]

    if use_noise:
        intrinsic_gal_ellip = 0.21  # The standard deviation of the intrinsic galaxy ellipticity distribution

        avg_gal_den = 30  # This is the average surface galaxy density in [num gals / arc min^2]
        area_per_pix = 1.49E8 / n_pix  # This is the total area in arcmin^2 divided by the number of pixels
        num_gal_per_pix = avg_gal_den * area_per_pix

        # Generate random Gaussian noise that will be added to our maps
        random_noise = np.random.normal(loc=0, scale=intrinsic_gal_ellip / np.sqrt(num_gal_per_pix), size=n_pix)

        # Compute what the expected Cl VALUE (singular) is for the shape noise.
        theory_cl_noise = intrinsic_gal_ellip ** 2 / (avg_gal_den / (sciconst.arcminute ** 2))

        # Include the l(l+1) / 2pi factor in the Cl values
        theory_cl_noise = ells * (ells + 1) * theory_cl_noise / (2 * np.pi)

    else:
        # If not using noise, then theory cl values are zero, along with the random noise elements
        theory_cl_noise = np.zeros(lmax - 1)

        random_noise = np.zeros(n_pix, dtype=np.bool)

    # Initiate a LCDM cosmology with a default cosmology
    params = camb.CAMBparams()
    params.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.06)
    params.InitPower.set_params(As=2.1E-9, ns=0.96)
    params.set_for_lmax(lmax, lens_potential_accuracy=0)

    params.Want_CMB = False
    params.WantCls = True
    params.NonLinear = camb.model.NonLinear_both

    # We want to evaluate the lensing power spectrum at z=2 only, for the moment.
    params.SourceWindows = [GaussianSourceWindow(redshift=0.5, source_type='lensing', sigma=0.05),
                            GaussianSourceWindow(redshift=2.0, source_type='lensing', sigma=0.05)]

    # Read in the unmasked convergence map that was created with A_s = 2.1E-9 and m_nu = 0.06
    converg_map1 = hp.read_map('KappaGammaMap-f2z1.fits', verbose=True, field=0)
    converg_map2 = hp.read_map('KappaGammaMap-f2z2.fits', verbose=True, field=0)

    # Read in the Euclid mask
    euclid_mask = hp.read_map('/home/maraio/Codes/LensingMapMaking/resources/Euclid_masks/Euclid-gal-mask-2048.fits',
                              verbose=False).astype(np.bool)

    f_sky = euclid_mask.sum() / euclid_mask.size

    # Add the random shape noise to our convergence map, this is zeros if not using noise
    converg_map1 += random_noise
    converg_map2 += random_noise

    # Create masked map from convergence map and Euclid mask
    converg_map_mask1 = hp.ma(converg_map1)
    converg_map_mask1.mask = np.logical_not(euclid_mask)

    converg_map_mask2 = hp.ma(converg_map2)
    converg_map_mask2.mask = np.logical_not(euclid_mask)

    # Convert to C_ells
    converg_cls1 = np.array(hp.anafast(converg_map1, lmax=lmax)[2:])
    converg_cls1 = ells * (ells + 1) * converg_cls1 / (2 * np.pi)

    converg_cls2 = np.array(hp.anafast(converg_map2, lmax=lmax)[2:])
    converg_cls2 = ells * (ells + 1) * converg_cls2 / (2 * np.pi)

    converg_cls_mask1 = np.array(hp.anafast(converg_map_mask1, lmax=lmax)[2:])
    converg_cls_mask1 = ells * (ells + 1) * converg_cls_mask1 / (2 * np.pi)

    converg_cls_mask2 = np.array(hp.anafast(converg_map_mask2, lmax=lmax)[2:])
    converg_cls_mask2 = ells * (ells + 1) * converg_cls_mask2 / (2 * np.pi)

    # Normalise masked Cl through 1/f_sky
    converg_cls_mask1 /= f_sky
    converg_cls_mask2 /= f_sky

    # Whether to use the masked Cls or not

    if options[option_section, "use_mask"]:
        # Masked sky
        loaded_data = {'params': params, 'converg_cls1': converg_cls_mask1, 'converg_cls2': converg_cls_mask2,
                       'ells': ells, 'theory_cl_noise': theory_cl_noise}
    else:
        # Unmasked sky
        loaded_data = {'params': params, 'converg_cls1': converg_cls1, 'converg_cls2': converg_cls2, 'ells': ells,
                       'theory_cl_noise': theory_cl_noise}

    # Return our loaded_data dictionary
    return loaded_data


def execute(block, config):
    loaded_data = config

    # Load in the ell values
    ells = loaded_data['ells']

    # Load in the CAMB parameter class
    params = loaded_data['params']

    # Load in the convergence cl values for the two redshift bins
    converg_cls1 = loaded_data['converg_cls1']
    converg_cls2 = loaded_data['converg_cls2']

    # Read in the theory cl noise, if we're not using noise this is just zeros
    theory_cl_noise = loaded_data['theory_cl_noise']

    # Read in the cosmological parameters for this specific run
    A_s = block[cosmo, "A_s"]
    n_s = block[cosmo, "n_s"]
    H0 = block[cosmo, "H0"]
    omega_b = block[cosmo, "omega_b"]
    omega_c = block[cosmo, "omega_c"]
    m_nu = block[cosmo, 'm_nu']
    w0 = block[cosmo, "w0"]
    wa = block[cosmo, "wa"]

    # Set the cosmological parameters in CAMB
    params.InitPower.set_params(As=A_s, ns=n_s)
    params.set_cosmology(H0=H0, ombh2=omega_b, omch2=omega_c, mnu=m_nu)
    params.DarkEnergy = DarkEnergyPPF(w=w0, wa=wa)

    # Compute the lensing power spectrum
    results = camb.get_results(params)

    # Get the lensing power spectrum Cl
    cl = results.get_source_cls_dict()

    # Separate the two different redshift bins
    cl_theory1 = cl['W1xW1'][2:2001]
    cl_theory2 = cl['W2xW2'][2:2001]

    # Compute the log-likelihoods for each redshift bin
    log_lik1 = -1 * np.sum(
        (2 * ells + 1) * (np.log(cl_theory1 + theory_cl_noise) + converg_cls1 / (cl_theory1 + theory_cl_noise)))
    log_lik2 = -1 * np.sum(
        (2 * ells + 1) * (np.log(cl_theory2 + theory_cl_noise) + converg_cls2 / (cl_theory2 + theory_cl_noise)))

    # Combine the two log-likelihoods
    log_lik = log_lik1 + log_lik2

    # Put the likelihood into the Cosmosis block
    block.put("likelihoods", "me_like", log_lik)

    # We tell Cosmosis that everything went fine by returning zero
    return 0


def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass

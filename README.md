Pulse Portraiture + Polarimetry
=================


## What?

A set of libraries and modules to measure "wideband" pulse times-of-arrival (TOAs), written in python. It uses an extension of Joe Taylor's **FFTFIT** algorithm (**Taylor 1992**) to simultaneously measure a phase (TOA) and dispersion measure (DM).  It has subsequently been improved to also incoporate fitting for scattering parameters (timescale tau and index alpha), frequency**-4 phase delays ("GM")**, and extened to work with all four Stokes parameters. It is to be used with [PSRCHIVE][psrchive]-compatible folded archives ([PSRFITS][psrfits] format).



## Why?

The motivation behind writing this software was to develop a wideband measurement routine for high-precision pulsar timing in the era of very broadband receivers and high-cadence timing observations for PTA experiments, and it makes up a chunk of Tim Pennucci's Ph.D. thesis. Algorithm development and coding help was provided by Paul Demorest and Scott Ransom.

## How?

The technical description of this work and its related papers are:

* [Pennucci, Demorest, & Ransom (2014), "_Elementary Wideband Timing of Radio Pulsars_", ApJ, 790, 93][2014].

  > **NB**: Equation **13** of the paper was written incorrectly: It should be the reverse:
  >
  > DM<sub>bary</sub> = DM<sub>topo</sub> / doppler_factor.

* [Pennucci (2015), "_Wideband Observations of Radio Pulsars_", PhDT, UVa][2015].
* [Pennucci (2019), "_Frequency-dependent Template Profiles for High-precision Pulsar Timing_", ApJ, 871, 1][2019].
* [Curylo (2025), "_Frequency- and phase-resolved polarimetry of millisecond pulsars and its application to timing_", arXiv:2512.09220][2025].

## Requirements

* [**PSRCHIVE**][psrchive], compiled with the python-interface enabled,
* [**NumPy**][numpy] & [**SciPy**][scipy] recent versions will do,
* [**PyWavelets**][pywt] is required for wavelet smoothing / ppspline.py, and
* [**LMFIT**][lmfit] is required for Gaussian portrait modelling / ppgauss.py and a few other functions.

## TL;DR

* [`pplib_pol`][pplib_pol] contains functions and classes needed for the fitting scripts.
* [`ppspline_pol`][ppspline_pol] is a command-line utility to build smoothly varying model portraits based on PCA decomposition, wavelet smoothing, and B-splin einterpolation between the components.
* [`ppgauss`][ppgauss] is a command-line utility to build Gaussian-component model portraits.
* [`pptoaslib`][pptoaslib] contains functions needed for pptoas.
* [`pptoas`][pptoas] is a command-line utility to measure TOAs, DMs, nu**-4 delays, and scattering parameters.
* [`ppalign_pol`][ppalign_pol] is a command-line utility to average homogeneous data by measuring phases and DMs.
* [`ppzap`][ppzap] is a command-line utility which uses pptoas to identify potentially overlooked bad channels to zap.
* The command-line programs can be imported into ipython for additional flexibility of use.
* See the [**examples**][examples] directory for simple command-line use.
* Run and examine [**examples/**`example.py`][examplepy] for a more in-depth demonstration (total intensity version)
* Try the jupyter notebooks for a walk-through in creating portraits, templates and polarimetry. 

## License

Released under **GPLv2**, sans "or later" clause.

## Other

Code improvements are underway, as is a broad application to IPTA pulsars of interest. [Suggestions and additional development are welcome](https://github.com/gcurylo/PulsePortraiture).

[psrfits]: https://www.atnf.csiro.au/research/pulsar/psrfits_definition/Psrfits.html

[2014]: https://doi.org/10.1088/0004-637X/790/2/93
[2015]: https://doi.org/10.18130/V3W56C
[2019]: https://doi.org/10.3847/1538-4357/aaf6ef
[2025]: https://arxiv.org/abs/2512.09220

[psrchive]: http://psrchive.sourceforge.net/
[numpy]: https://numpy.org/
[scipy]: https://www.scipy.org/
[pywt]: https://pywavelets.readthedocs.io/en/latest/ref/wavelets.html
[lmfit]: https://lmfit.github.io/lmfit-py/index.html

[pplib_pol]: https://github.com/gcurylo/PulsePortraiture/blob/master/pplib_pol.py
[ppspline_pol]: https://github.com/gcurylo/PulsePortraiture/blob/master/ppspline_pol.py
[ppgauss]: https://github.com/gcurylo/PulsePortraiture/blob/master/ppgauss.py
[pptoaslib]: https://github.com/gcurylo/PulsePortraiture/blob/master/pptoaslib.py
[pptoas]: https://github.com/gcurylo/PulsePortraiture/blob/master/pptoas.py
[ppalign_pol]: https://github.com/gcurylo/PulsePortraiture/blob/master/ppalign_pol.py
[ppzap]: https://github.com/gcurylo/PulsePortraiture/blob/master/ppzap.py
[examples]: https://github.com/gcurylo/PulsePortraiture/tree/master/examples
[examplepy]: https://github.com/gcurylo/PulsePortraiture/blob/master/examples/example.py

#!/usr/bin/env python

from setuptools import setup

setup(name='PulsePortraiture',
      version='0.2.0',
      description='Data analysis package for wideband pulsar timing including polarimetry.',
      author='Tim Pennucci (main), Gosia Curylo (polarimetry)',
      author_email='tim.pennucci@nanograv.org, gosia.curylo@monash.edu',
      url='https://github.com/gcurylo/PulsePortraiture',
      py_modules=['ppalign_pol', 'ppgauss', 'pplib_pol', 'ppspline_pol', 'pptoas', 'pptoaslib', 'ppzap','telescope_codes'],
      scripts=['ppalign_pol.py','ppgauss.py', 'ppspline_pol.py', 'pptoas.py',
          'ppzap.py']
     )

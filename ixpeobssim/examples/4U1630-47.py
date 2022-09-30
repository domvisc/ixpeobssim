from __future__ import print_function, division
from tkinter import N

import os

from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.config.microquasar_4U1630 as input_model
from ixpeobssim import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_



DURATION = 460000.
expression = 'tbabs*(diskbb+powerlaw)' #xspec spectrum model 
paramsfile = os.path.join(IXPEOBSSIM_CONFIG_ASCII,'4U1630-47_fit_params.csv')


def simulate():
    """Run the simulation and fold the events in phase.
    """
    pipeline.xpobssim(duration=DURATION,configfile='/home/dom/Desk/Tesi/ixpeobssim/ixpeobssim/config/microquasar_4U1630.py',overwrite =True)


def bin_():
    """Create the pha1 files.
    """
    for algorithm in ['PCUBE']:
        pipeline.xpbin(*pipeline.file_list(), algorithm=algorithm, irfname = 'ixpe:obssim:v11', overwrite = True, ebinalg = 'LIN' , ebins = 4)


def spectro_polarimetric_fit(expression,paramsfile):
    """Perform a two-tier spectro-polarimetric fit in XSPEC.
    """
    if not PYXSPEC_INSTALLED:
        return

    file_list = pipeline.file_list('pha1*')
    model = '%s * pollin' % expression
    fit_output = pipeline.xpxspec(*file_list, model=expression, paramsfile=paramsfile, plot=True, error=False)


def draw_spectral_model(expression,parameters):

    parameters = xspec_.read_parameter_file(paramsfile)
    binning, energy, flux, parameter_names = xspec_.sample_spectral_model(expression,parameters,emin=1.,emax=12.,num_points=1000)
    plt.figure('input model')
    plt.plot(energy,flux,'o')
    


def run():
    """Run all.
    """
    simulate()
    bin_()
    #draw_spectral_model(expression,paramsfile)
    #spectro_polarimetric_fit(expression,paramsfile)



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('microquasar_4U1630')

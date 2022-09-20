from __future__ import print_function, division
from tkinter import N

import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.config.microquasar_4U1630 as input_model
from ixpeobssim import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_


DURATION = 50000.

def simulate():
    """Run the simulation and fold the events in phase.
    """
    pipeline.xpobssim(duration=DURATION)


def bin_():
    """Create the pha1 files.
    """
    for algorithm in ['PHA1', 'PHA1Q', 'PHA1U']:
        pipeline.xpbin(*pipeline.file_list(), algorithm=algorithm)


def spectro_polarimetric_fit():
    """Perform a two-tier spectro-polarimetric fit in XSPEC.
    """
    if not PYXSPEC_INSTALLED:
        return
    
    expression = 'phabs * (bbody + powerlaw)'

    col_dens = 7.8 #e22
    bb_temp = 1.38 #kev
    bb_norm = 158
    pl_index = 1.7
    pl_norm = 1

    file_list = pipeline.file_list('pha1*')
    model = '%s * polconst' % expression
    fit_output = pipeline.xpxspec(*file_list, model=model, plot=True)
    target = (col_dens, bb_temp, bb_norm,
              pl_index, pl_norm, input_model.pol_deg,
              input_model.pol_ang)
    xspec_.compare_fit_data(fit_output, target)



def run():
    """Run all.
    """
    bin_()
    spectro_polarimetric_fit()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('microquasar_4U1630')

import os

import numpy

from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
from ixpeobssim.core.spline import xInterpolatedUnivariateLogSpline
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import ergcms_to_mcrab
from ixpeobssim.srcmodel.polarization import constant

__model__ = file_path_to_model_name(__file__)

SRC_NAME = '4U1630-47'
SRC_RA = 248.506708
SRC_DEC = -47.393

def pol_deg(E,t=None,ra=None,dec=None):
    """Simple linear polarization model from quick look data"""
    return 0.044 + 0.0076*E

pol_ang = constant(numpy.radians(60))


def _parse_spec(file_name):
    """Parse the spectral data.
    Digitized data point from NUSTAR observation id 80802313002, 25/08/2022, 
    reduced by Andrea Marinucci
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, file_name)
    logger.info('Loading spectral data from %s...', file_path)
    energy, flux = numpy.loadtxt(file_path, unpack=True, delimiter=',')
    spec_spline = xInterpolatedUnivariateLogSpline(energy, flux, k=1, **fmtaxis.spec)
    return spec_spline

spec_spline = _parse_spec('4U1630-47_spectrum.csv')
spec = lambda E, t=None: spec_spline(E)
source = xPointSource(SRC_NAME,SRC_RA,SRC_DEC,spec,pol_deg,pol_ang)

ROI_MODEL = xROIModel(SRC_RA,SRC_DEC,source)

def display_spectrum(energy):
    """Plot the energy spectrum of the source.
    """
    # Calculate the standard integral flux.
    emin = 2.
    emax = 10.
    source_dict = {name: src for name, src in ROI_MODEL.items() if src.identifier < 3}
    flux_ergcms = sum(src.calculate_integral_flux(emin, emax, t=0) for src in source_dict.values())
    flux_mcrab = ergcms_to_mcrab(flux_ergcms)
    flux_label = 'Total flux: %.2e erg cm$^{-2}$ s$^{-1}$ (%.1f mcrab)' % (flux_ergcms, flux_mcrab)

    plt.figure('%s spectrum' % __model__)
    for name, component in source_dict.items():
        plt.plot(energy, component.photon_spectrum(energy), label=name)
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=1.e-3, logx=True,
              logy=True, grids=True, legend=True, **fmtaxis.spec)

def display_polarization(energy):
    """Plot the polarization degree and angle as a function of the energy.
    """
    source_dict = {name: src for name, src in ROI_MODEL.items() if src.identifier < 3}
    plt.figure('%s polarization degree' % __model__)
    for name, component in source_dict.items():
        plt.plot(energy, component.polarization_degree(energy), label=name)
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=-0.001, ymax=0.15,
              legend=True, grids=True, **fmtaxis.ene_pol_deg)

    plt.figure('%s polarization angle' % __model__)
    for name, component in source_dict.items():
        plt.plot(energy, numpy.degrees(component.polarization_angle(energy)), label=name)
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=-90., ymax=90.,
              legend=True, grids=True, **fmtaxis.ene_pol_ang_deg)


def display(emin=1., emax=15.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 1000)
    display_spectrum(energy)
    display_polarization(energy)


if __name__ == '__main__':
    bootstrap_display()
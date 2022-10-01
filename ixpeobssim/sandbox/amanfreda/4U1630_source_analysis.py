from __future__ import print_function, division
from tkinter import N

import numpy
import os
from functools import reduce

from ixpeobssim.binning.misc import xBinnedLightCurve
from ixpeobssim.utils.matplotlib_ import plt, color_wheel_ar, color_wheel_mpr
from ixpeobssim.utils.time_ import met_to_mjd
from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.evt.xspec_ as xspec_
from ixpeobssim.core.hist import xScatterPlot


DATA_FOLDER = '/home/alberto/xpe/xpedata/4U1630-472/01250401/event_l2'
file_names = ['ixpe01250401_det1_evt2_v03.fits',
              'ixpe01250401_det2_evt2_v03.fits',
              'ixpe01250401_det3_evt2_v03.fits']
file_list = [os.path.join(DATA_FOLDER, file_name) for file_name in file_names]


def combine_light_curves(lc1, lc2):
    """
    """
    lc1 += lc2
    return lc1

def apply_exposure_mask(lc, min_exposure=60):
    """
    """
    mask = lc.EXPOSURE < min_exposure
    lc.COUNTS[mask] = 0
    lc.ERROR[mask] = 0.
    lc.EXPOSURE[mask] = 0.

def draw_light_curve(light_curve, suffix='LC', label=None):
    plt.figure(suffix)
    light_curve.plot(mjd=True, label=label)
    if label is not None:
        plt.legend()
    plt.grid(True)

def draw_hardness_ratio(lowe_lc, highe_lc, label=None, mjd=True):
    """
    """
    rate_lowe = lowe_lc.rate()
    rate_highe = highe_lc.rate()
    rate_err_lowe = lowe_lc.rate_error()
    rate_err_highe = highe_lc.rate_error()
    plt.figure('hardness_ratio')
    ratio = (rate_highe - rate_lowe)/ (rate_highe + rate_lowe)
    ratio_err = 2 * numpy.sqrt(rate_lowe**2 * rate_err_highe**2 +
                               rate_highe**2 * rate_err_lowe**2) \
                / (rate_lowe + rate_highe)**2
    if mjd:
        t = met_to_mjd(lowe_lc.TIME)
        xlabel = 'MJD'
    else:
        t = lowe_lc.TIME
        xlabel = 'MET [s]'
    xScatterPlot(t, ratio, ratio_err, xlabel=xlabel,
                 ylabel='hardness ratio').plot(label=label)
    plt.grid(True)
    if label is not None:
        plt.legend()

def draw_time_regions(*boundaries, mjd=True):
    if mjd:
        boundaries = [met_to_mjd(t) for t in boundaries]
    for i, (t0, t1) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        plt.gca().axvspan(t0, t1, alpha=0.2, color=color_wheel_mpr(i))


def build_light_curve(suffix='eall', emin=2., emax=8., tbins=230):
    eselect_files = pipeline.xpselect(*file_list, suffix=suffix, emin=emin,
                                      emax=emax)
    lc_files = pipeline.xpbin(*eselect_files, algorithm='LC', tbins=tbins,
                              irfname='v11:ixpe:obssim')
    light_curves = [xBinnedLightCurve(lc_file) for lc_file in lc_files]
    light_curve = reduce(combine_light_curves, light_curves)
    apply_exposure_mask(light_curve)
    draw_light_curve(light_curve, suffix=suffix, label=f'[{emin:.1f}--{emax:.1f}] keV')
    return light_curve

def select_time_regions(boundaries, suffix='eall')

def run():
    """
    """
    tbins = 100
    all_energy_lc = build_light_curve(suffix='eall', emin=2., emax=8., tbins=tbins)
    low_energy_lc = build_light_curve(suffix='elow', emin=2., emax=3., tbins=tbins)
    high_energy_lc = build_light_curve(suffix='ehigh', emin=5., emax=8., tbins=tbins)
    draw_hardness_ratio(low_energy_lc, high_energy_lc,
                        label='[5-8] keV / [2-3] keV', mjd=True)
    t = all_energy_lc.TIME
    boundaries = [min(t), 178288000., 178347500., 178742060., max(t)]
    draw_time_zone(boundaries, mjd=True)
    select_time_regions(boundaries)



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('4U1630')

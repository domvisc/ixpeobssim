from __future__ import print_function, division
from tkinter import N

import numpy
import os
from functools import reduce

from ixpeobssim.binning.misc import xBinnedLightCurve
from ixpeobssim.utils.matplotlib_ import plt, color_wheel_ar, color_wheel_mpr
from ixpeobssim.utils.time_ import met_to_mjd
import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.evt.xspec_ as xspec_
from ixpeobssim.core.hist import xScatterPlot


DATA_FOLDER = '/home/alberto/xpe/xpedata/4U1630-472/01250401/event_l2'
FILE_NAMES = ['ixpe01250401_det1_evt2_v03.fits',
              'ixpe01250401_det2_evt2_v03.fits',
              'ixpe01250401_det3_evt2_v03.fits']
FILE_LIST = [os.path.join(DATA_FOLDER, file_name) for file_name in FILE_NAMES]


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

def draw_light_curve(light_curve, suffix=None, label=None):
    if suffix is None:
        suffix = 'LC'
    else:
        suffix = f'LC_{suffix}'
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

def draw_time_boundaries(boundaries, mjd=True):
    if mjd:
        boundaries = [met_to_mjd(t) for t in boundaries]
    for i, (t0, t1) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        plt.gca().axvspan(t0, t1, alpha=0.2, color=color_wheel_mpr(i))

def energy_select(*file_list, emin=2., emax=8., suffix=None):
    return pipeline.xpselect(*file_list, suffix=suffix, emin=emin, emax=emax)

def time_select(*file_list, tmin=None, tmax=None, suffix=None):
    if tmin is None:
        tmin = 0.
    if tmax is None:
        tmax = 100000000000.
    return pipeline.xpselect(*file_list, suffix=suffix, tmin=tmin, tmax=tmax)

def build_light_curve(*file_list, tbins=230, suffix=None, label=None):
    lc_files = pipeline.xpbin(*file_list, algorithm='LC', tbins=tbins,
                              irfname='v11:ixpe:obssim')
    light_curves = [xBinnedLightCurve(lc_file) for lc_file in lc_files]
    light_curve = reduce(combine_light_curves, light_curves)
    apply_exposure_mask(light_curve)
    draw_light_curve(light_curve, suffix=suffix, label=label)
    return light_curve

def time_suffix(i):
    return f'timesel_{i}'

def elabel(emin, emax):
    return f'[{emin:.1f}--{emax:.1f}] keV'

def select_time_regions(*file_list, boundaries=None):
    if boundaries is None:
        return
    for i, (t0, t1) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        pipeline.xpselect(*file_list, tmin=t0, tmax=t1, suffix=time_suffix(i))

def run():
    """
    """
    tbins = 100
    all_energy_files = energy_select(*FILE_LIST, suffix='eall', emin=2., emax=8.)
    low_energy_files = energy_select(*FILE_LIST, suffix='elow', emin=2., emax=3.)
    hig_energy_files = energy_select(*FILE_LIST, suffix='ehigh', emin=5., emax=8.)
    all_energy_lc = build_light_curve(*all_energy_files, suffix='eall',
                                      label=elabel(2., 8.), tbins=tbins)
    t = all_energy_lc.TIME
    boundaries = [min(t), 178288000., 178347500., 178742060., max(t)]
    low_energy_lc = build_light_curve(*low_energy_files, suffix='elow',
                                      label=elabel(2., 3.),tbins=tbins)
    high_energy_lc = build_light_curve(*hig_energy_files, suffix='ehigh',
                                      label=elabel(5., 8.),tbins=tbins)
    draw_hardness_ratio(low_energy_lc, high_energy_lc,
                        label='[5-8] keV / [2-3] keV', mjd=True)
    t = all_energy_lc.TIME
    boundaries = [min(t), 178288000., 178347500., 178742060., max(t)]
    draw_time_boundaries(boundaries, mjd=True)
    print('Selecting time zones...')
    print(all_energy_files)
    input()
    select_time_regions(*all_energy_files, boundaries=boundaries)



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('4U1630')

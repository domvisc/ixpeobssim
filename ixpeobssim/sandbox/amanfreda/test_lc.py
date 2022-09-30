import numpy
import os

from ixpeobssim.binning.misc import xBinnedLightCurve
from ixpeobssim.utils.matplotlib_ import plt

def find_saa_edges(saa_mets):
    saa_switch = numpy.ediff1d(saa_mets) < 10.
    saa_edges_mask = numpy.full(len(saa_mets), False, dtype=bool)
    saa_edges_mask[0] = True
    saa_edges_mask[-1] = True
    saa_edges_mask[1:] += saa_switch
    saa_edges_mask[:-1] += saa_switch
    saa_edges = saa_mets[saa_edges_mask]
    return saa_edges

def apply_saa_mask(lc, saa_edges, saa_dist=60):
    """
    """
    idx = numpy.searchsorted(saa_edges, lc.TIME)
    mask = numpy.full(lc.TIME.shape, False, dtype=bool)
    for i, saa_idx in enumerate(idx):
        t = lc.TIME[i]
        if saa_idx == 0:
            mask[i] = numpy.abs(saa_edges[saa_idx] - t) < saa_dist
        if saa_idx == len(saa_edges):
            mask[i] = numpy.abs(saa_edges[saa_idx - 1] - t) < saa_dist
        else:
            mask[i] = (numpy.abs(saa_edges[saa_idx - 1] - t) < saa_dist) | \
                      (numpy.abs(saa_edges[saa_idx] - t) < saa_dist)
    lc.COUNTS[mask] = 0
    lc.ERROR[mask] = 0.
    lc.EXPOSURE[mask] = 0.

def apply_exposure_mask(lc, min_exposure=30):
    """
    """
    mask = lc.EXPOSURE < min_exposure
    lc.COUNTS[mask] = 0
    lc.ERROR[mask] = 0.
    lc.EXPOSURE[mask] = 0.


DATA_FOLDER = '/home/alberto/xpe/xpedata/4U1630-472/01250401/event_l2'
saa_mets = numpy.load('/home/alberto/xpe/xpedata/4U1630-472/4U_1630_saa_mets.npz')['saa_mets']
#saa_edges = find_saa_edges(saa_mets)

def plot_lc_separate(*lcs, label='[2-8] keV', title='LC', draw_saa=False):
    """
    """
    for lc in lcs:
        apply_exposure_mask(lc)
    fig = plt.figure(f'{title}_separate', figsize=(10,6))
    for lc in lcs:
        lc.plot(mjd=True, label=f'DU{lc.du_id()} {label}')
    ymin, ymax = plt.gca().get_ylim()
    plt.ylim(0., ymax)
    plt.grid(True)
    if draw_saa:
        ax = plt.gca()
        ax.fill_between(saa_mets[:-1], 0, 25, where=numpy.ediff1d(saa_mets) < 10.,
                        color='gray', alpha=0.5, transform=ax.get_xaxis_transform(),
                        label='SAA')
    plt.legend()

def plot_lc(*lcs, label='[2-8] keV', title='LC', draw_saa=False):
    """
    """
    fig = plt.figure(title, figsize=(10,6))
    for lc in lcs[1:]:
        lcs[0] += lc
    lcs[0].plot(mjd=True, label=f'All DUs {label}')
    ymin, ymax = plt.gca().get_ylim()
    plt.ylim(0., ymax)
    plt.grid(True)
    if draw_saa:
        ax = plt.gca()
        ax.fill_between(saa_mets[:-1], 0, 25, where=numpy.ediff1d(saa_mets) < 10.,
                        color='gray', alpha=0.5, transform=ax.get_xaxis_transform(),
                        label='SAA')
    plt.legend()


alle_names = ['ixpe01250401_det1_evt2_v03_eselect_lc.fits',
              'ixpe01250401_det2_evt2_v03_eselect_lc.fits',
              'ixpe01250401_det3_evt2_v03_eselect_lc.fits']
highe_names = ['ixpe01250401_det1_evt2_v03_ehigh_lc.fits',
               'ixpe01250401_det2_evt2_v03_ehigh_lc.fits',
               'ixpe01250401_det3_evt2_v03_ehigh_lc.fits']
lowe_names = ['ixpe01250401_det1_evt2_v03_elow_lc.fits',
              'ixpe01250401_det2_evt2_v03_elow_lc.fits',
              'ixpe01250401_det3_evt2_v03_elow_lc.fits']

alle_lcs = [xBinnedLightCurve.from_file_list([os.path.join(DATA_FOLDER,
            file_name)]) for file_name in alle_names]
lowe_lcs = [xBinnedLightCurve.from_file_list([os.path.join(DATA_FOLDER,
            file_name)]) for file_name in lowe_names]
highe_lcs = [xBinnedLightCurve.from_file_list([os.path.join(DATA_FOLDER,
             file_name)]) for file_name in highe_names]

for lc in alle_lcs[1:]:
    alle_lcs[0] += lc

for lc in lowe_lcs[1:]:
    lowe_lcs[0] += lc

for lc in highe_lcs[1:]:
    highe_lcs[0] += lc

apply_exposure_mask(alle_lcs[0], min_exposure=60)
apply_exposure_mask(highe_lcs[0], min_exposure=60)
apply_exposure_mask(lowe_lcs[0], min_exposure=60)

plt.figure('LC')
alle_lcs[0].plot(mjd=True)


rate_lowe = lowe_lcs[0].rate()
rate_highe = highe_lcs[0].rate()
rate_err_lowe = lowe_lcs[0].rate_error()
rate_err_highe = highe_lcs[0].rate_error()

plt.figure('ratio')
ratio = rate_lowe/rate_highe
ratio_err = numpy.sqrt(rate_err_lowe**2 + rate_err_highe**2 * rate_lowe) / rate_highe
plt.errorbar(lowe_lcs[0].TIME, ratio, yerr=ratio_err, fmt='o',
             label='[2-3] keV / [5-8] keV')
plt.grid(True)
plt.xlabel('MET [s]')
plt.ylabel('rate ratio')
plt.legend()
#plot_lc(*highe_names, label='[5-8] keV', title=None, draw_saa=False)
#plot_lc(*lowe_names, label='[2-3] keV', title=None, draw_saa=False)

plt.ion()
plt.show()

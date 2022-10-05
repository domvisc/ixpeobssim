from __future__ import print_function, division
from cProfile import label
from tkinter import N

import numpy
import os
import logging

from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
from ixpeobssim.evt.event import xEventFile
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.polarization import xBinnedPolarizationCube


input_files=['ixpe01250401_det1_evt2_v03.fits',
             'ixpe01250401_det2_evt2_v03.fits',
             'ixpe01250401_det3_evt2_v03.fits']
RA = 248.5042
DEC = -47.3900


def set_work_directory_and_filter_in_energy(is_the_path_set = False ,is_filtered_in_energy = False):
    '''Set the directory of the .fits file, if is_filtered_in_energy is False they are filtered
    (and overwrited) in energy with xpselect ---- 2 - 8 keV.
    '''
    if is_the_path_set is False:
        work_dir = input('write path (str) to input_files dir, they will be overwrited if is_filtered_in_energy is False! \n')
        os.chdir(work_dir)
    else:
         os.chdir('/home/dom/Desk/Tesi/Data/IXPE_DATA/') #setting work directory manually

    if is_filtered_in_energy is False:
        pipeline.xpselect(*input_files, emin = 2, emax = 8)

WORK_DIR = '/home/alberto/xpe/xpedata/4U1630-472/01250401/event_l2'
input_files = [os.path.join(WORK_DIR, f.replace('.fits', '_eall.fits')) \
               for f in input_files]

def create_data_file():
    ''' create data file at distance from the approximate center of [1,1.5,2,2.5,3.] arcmin and cmaps
    '''
    for RAD in ([1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,6.,7.]):

        bkg_file = pipeline.xpselect(*input_files, suffix = f'bkg_{RAD}', innerrad = RAD, ra = RA, dec = DEC)
        data_file = pipeline.xpselect(*input_files, suffix = f'data_{RAD}', rad = RAD, ra = RA, dec = DEC)
        pipeline.xpbin(*data_file, algorithm = 'CMAP', irfname = 'ixpe:obssim:v11')
        pipeline.xpbin(*bkg_file, algorithm = 'CMAP', irfname = 'ixpe:obssim:v11')
        pipeline.xpbin(*data_file, algorithm = 'PCUBE', irfname = 'ixpe:obssim:v11', ebins = 4)
        pipeline.xpbin(*bkg_file, algorithm = 'PCUBE', irfname = 'ixpe:obssim:v11', ebins = 4)
        pipeline.xpbin(*data_file, algorithm = 'PMAPCUBE', irfname = 'ixpe:obssim:v11', ebins = 1)
        pipeline.xpbin(*bkg_file, algorithm = 'PMAPCUBE', irfname = 'ixpe:obssim:v11', ebins = 1)

def data_files(rad, suffix=None):
    if suffix is None:
        suffix = f'{rad}'
    else:
        suffix = f'{rad}_{suffix}'
    return [f.replace('.fits', f'_data_{suffix}.fits') for f in input_files]

def bkg_files(rad, suffix=None):
    if suffix is None:
        suffix = f'{rad}'
    else:
        suffix = f'{rad}_{suffix}'
    return [f.replace('.fits', f'_bkg_{suffix}.fits') for f in input_files]

def evt_fraction(rad_list):
    data_evts = []
    bkg_evts = []
    evt_fraction_dict = {}
    for RAD in rad_list:
        _data_evts = 0
        _bkg_evts = 0
        data_paths = data_files(RAD)
        bkg_paths = bkg_files(RAD)
        for data_path in data_paths:
            evt_file = xEventFile(data_path)
            _data_evts += evt_file.num_events()
        data_evts.append(_data_evts)
        for bkg_path in bkg_paths:
            evt_file = xEventFile(bkg_path)
            _bkg_evts += evt_file.num_events()
        bkg_evts.append(_bkg_evts)
        evt_fraction_dict[RAD] = _data_evts / (_data_evts + _bkg_evts)
    plt.figure('evt_fraction_vs_radius')
    plt.plot(rad_list, evt_fraction_dict.values(), 'o')
    plt.xlabel('Enclosing radius [arcmin]')
    plt.ylabel('Evt. fraction')
    plt.grid(True)
    return evt_fraction_dict

def view_cmaps(rad_list):
    '''visualize cmaps given a rad list'''

    for RAD in rad_list:

        data_path = data_files(RAD, 'CMAP')
        bkg_path = bkg_files(RAD, 'CMAP')

        pipeline.xpbinview(*data_path,stretch = 'log') #sum of dus?
        pipeline.xpbinview(*bkg_path,stretch = 'log')

def view_pcubes(rad_list):
    '''visualize pcubes given a rad list'''
    for RAD in rad_list:

        data_path = data_files(RAD, 'pcube')
        bkg_path = bkg_files(RAD, 'pcube')

        pipeline.xpbinview(*data_path) #sum of dus?
        pipeline.xpbinview(*bkg_path)

def plot_pol_vs_energy(energy, pd, pd_err, label=None, title=None, **ax_setup):
    plt.errorbar(energy, pd, yerr=pd_err, fmt='o', label=label, linestyle='--')
    if label is not None:
        legend = True
    else:
        legend = False
    if title is not None:
        plt.title(title, fontsize='large', fontweight='semibold')
    setup_gca(legend=legend, **ax_setup)

def read_pcubes(rad_list, evt_fract_dict=None):
    pol_deg_plot_setup = {'xlabel' : 'E [keV]',
                          'ylabel' : 'polarization degree',
                          'grids'  : True,
                          'ymin'   : 0.03,
                          'ymax'   : 0.15}
    pol_ang_plot_setup = {'xlabel' : 'E [keV]',
                          'ylabel' : 'polarization angle [deg]',
                          'grids'  : True,
                          'ymin'   : 12.,
                          'ymax'   : 24.}
    for RAD in rad_list:
        data_path = data_files(RAD, 'pcube')
        bkg_path = bkg_files(RAD, 'pcube')
        label = f'radius = {RAD} [arcmin]'
        if evt_fract_dict is not None:
            label = f'{label} ({100*evt_fract_dict[RAD]:.1f} % evts)'
        polarization_data_cube = xBinnedPolarizationCube.from_file_list(data_path)
        _pd, _pd_err, _pa, _pa_err = polarization_data_cube.polarization()
        plt.figure('pd_vs_rad')
        plot_pol_vs_energy(polarization_data_cube.E_MEAN, _pd, _pd_err,
                           label=label, title='NO BKG SUB', **pol_deg_plot_setup)
        plt.figure('pa_vs_rad')
        plot_pol_vs_energy(polarization_data_cube.E_MEAN, _pa, _pa_err,
                           label=label, title='NO BKG SUB', **pol_ang_plot_setup)
        polarization_bkg_cube = xBinnedPolarizationCube.from_file_list(bkg_path)
        _pd_bkg, _pd_err_bkg, _pa_bkg, _pa_err_bkg = polarization_bkg_cube.polarization()
        plt.figure('bkg_pd_vs_rad')
        plot_pol_vs_energy(polarization_data_cube.E_MEAN, _pd_bkg, _pd_err_bkg,
                           label=label, **pol_deg_plot_setup)
        plt.figure('bkg_pa_vs_rad')
        plot_pol_vs_energy(polarization_data_cube.E_MEAN, _pa_bkg, _pa_err_bkg,
                           label=label, **pol_ang_plot_setup)
        polarization_bkg_cube *= polarization_data_cube.backscal() / polarization_bkg_cube.backscal()
        polarization_data_cube -= polarization_bkg_cube
        _pd, _pd_err, _pa, _pa_err = polarization_data_cube.polarization()
        plt.figure('bkg_sub_pd_vs_rad')
        plot_pol_vs_energy(polarization_data_cube.E_MEAN, _pd, _pd_err,
                           label=label, title='BKG SUB', **pol_deg_plot_setup)
        plt.figure('bkg_sub_pa_vs_rad')
        plot_pol_vs_energy(polarization_data_cube.E_MEAN, _pa, _pa_err,
                          label=label, title='BKG SUB', **pol_ang_plot_setup)


def view_pmapcubes(rad_list):
    '''visualize pcubes given a rad list'''
    for RAD in rad_list:
        data_path = data_files(RAD, 'pmapcube')
        bkg_path = bkg_files(RAD, 'pmapcube')
        pipeline.xpbinview(*data_path) #sum of dus?
        pipeline.xpbinview(*bkg_path)

def run():
    #set_work_directory_and_filter_in_energy(is_the_path_set = True, is_filtered_in_energy = True)
    #create_data_file()
    #view_cmaps([1.,7.])
    #view_pcubes([1.,1.5,2.,2.5,3.,3.5])
    #view_pmapcubes([1., 2., 3.5])
    evt_fract_dict = evt_fraction([1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7.])
    read_pcubes([1.,2., 3., 4.], evt_fract_dict)



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('microquasar_4U1630')
    plt.show()

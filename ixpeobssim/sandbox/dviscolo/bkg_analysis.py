from __future__ import print_function, division
from cProfile import label
from tkinter import N

import os
import logging

from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.polarization import xBinnedPolarizationCube


input_files=['det_1.fits','det_2.fits','det_3.fits']
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
        #pipeline.xpbin(*data_file, algorithm = 'PMAPCUBE', irfname = 'ixpe:obssim:v11', ebins = 1)
        #pipeline.xpbin(*bkg_file, algorithm = 'PMAPCUBE', irfname = 'ixpe:obssim:v11', ebins = 1)


    
def view_cmaps(rad_list):
    '''visualize cmaps given a rad list'''

    for RAD in rad_list:
        
        data_path = [f'det_1_data_{RAD}_cmap.fits',f'det_2_data_{RAD}_cmap.fits',f'det_2_data_{RAD}_cmap.fits']
        bkg_path = [f'det_1_bkg_{RAD}_cmap.fits',f'det_2_bkg_{RAD}_cmap.fits',f'det_1_bkg_{RAD}_cmap.fits']
        
        pipeline.xpbinview(*data_path,stretch = 'log') #sum of dus?
        pipeline.xpbinview(*bkg_path,stretch = 'log')

def view_pcubes(rad_list):
    '''visualize pcubes given a rad list'''

    for RAD in rad_list:
        
        data_path = [f'det_1_data_{RAD}_pcube.fits',f'det_2_data_{RAD}_pcube.fits',f'det_2_data_{RAD}_pcube.fits']
        bkg_path = [f'det_1_bkg_{RAD}_pcube.fits',f'det_2_bkg_{RAD}_pcube.fits',f'det_1_bkg_{RAD}_pcube.fits']
        
        pipeline.xpbinview(*data_path) #sum of dus?
        pipeline.xpbinview(*bkg_path)

def read_pcubes(rad_list):

    pd = []
    pd_err = []


    for RAD in rad_list:
        
        data_path = [f'det_1_data_{RAD}_pcube.fits',f'det_2_data_{RAD}_pcube.fits',f'det_2_data_{RAD}_pcube.fits']
        bkg_path = [f'det_1_bkg_{RAD}_pcube.fits',f'det_2_bkg_{RAD}_pcube.fits',f'det_1_bkg_{RAD}_pcube.fits']
        polarization_data_cube = xBinnedPolarizationCube.from_file_list(data_path)
        _pd, _pd_err, _pa, _pa_err = polarization_data_cube.polarization()
        plt.figure('pd_vs_rad')
        plt.errorbar(polarization_data_cube.E_MEAN,_pd,yerr=_pd_err,fmt='o',label=f'RAD={RAD}',linestyle = '--')
        plt.xlabel('E [keV]')
        plt.ylabel('polarization degree')
        plt.legend()
        plt.grid(True)
        plt.figure('pa_vs_rad')
        plt.errorbar(polarization_data_cube.E_MEAN,_pa,yerr=_pa_err,fmt='o',label=f'RAD={RAD}',linestyle = '--')
        plt.xlabel('E [keV]')
        plt.ylabel('polarization angle [deg]')
        plt.legend()
        plt.grid(True)
        polarization_bkg_cube = xBinnedPolarizationCube.from_file_list(bkg_path)
        _pd_bkg, _pd_err_bkg, _pa_bkg, _pa_err_bkg = polarization_bkg_cube.polarization()
        plt.figure('bkg_pd_vs_rad')
        plt.errorbar(polarization_bkg_cube.E_MEAN,_pd_bkg,yerr=_pd_err_bkg,fmt='o',label=f'RAD={RAD}',linestyle = '--')
        plt.xlabel('E [keV]')
        plt.ylabel('polarization degree')
        plt.legend()
        plt.grid(True)
        plt.figure('bkg_pa_vs_rad')
        plt.errorbar(polarization_bkg_cube.E_MEAN,_pa_bkg,yerr=_pa_err_bkg,fmt='o',label=f'RAD={RAD}',linestyle = '--')
        plt.xlabel('E [keV]')
        plt.ylabel('polarization angle [deg]')
        plt.legend()
        plt.grid(True)
        polarization_data_cube -= polarization_bkg_cube
        _pd, _pd_err, _pa, _pa_err = polarization_data_cube.polarization()
        plt.figure('bkg_sub_pd_vs_rad')
        plt.errorbar(polarization_data_cube.E_MEAN,_pd,yerr=_pd_err,fmt='o',label=f'RAD={RAD}',linestyle = '--')
        plt.xlabel('E [keV]')
        plt.ylabel('polarization degree')
        plt.legend()
        plt.grid(True)
        plt.figure('bkg_pa_vs_rad')
        plt.figure('bkg_sub_pa_vs_rad')
        plt.errorbar(polarization_data_cube.E_MEAN,_pa,yerr=_pa_err,fmt='o',label=f'RAD={RAD}',linestyle = '--')
        plt.xlabel('E [keV]')
        plt.ylabel('polarization angle [deg]')
        plt.legend()
        plt.grid(True)
    

def view_pmapcubes(rad_list):
    '''visualize pcubes given a rad list'''

    for RAD in rad_list:
        
        data_path = [f'det_1_data_{RAD}_pmapcube.fits',f'det_2_data_{RAD}_pmapcube.fits',f'det_2_data_{RAD}_pmapcube.fits']
        bkg_path = [f'det_1_bkg_{RAD}_pmapcube.fits',f'det_2_bkg_{RAD}_pmapcube.fits',f'det_1_bkg_{RAD}_pmapcube.fits']
        
        pipeline.xpbinview(*data_path) #sum of dus?
        pipeline.xpbinview(*bkg_path)


def run():

    set_work_directory_and_filter_in_energy(is_the_path_set = True, is_filtered_in_energy = True)
    create_data_file()
    #view_cmaps([1.,7.])
    #view_pcubes([1.,1.5,2.,2.5,3.,3.5])
    #view_pmapcubes([1.,1.5,2.,2.5,3.,3.5])
    read_pcubes([1.,4.,7.])



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('microquasar_4U1630')
    plt.show()



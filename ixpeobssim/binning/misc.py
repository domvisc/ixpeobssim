# Copyright (C) 2015--2022, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""Miscellanea binned products.
"""

from __future__ import print_function, division

from astropy.io import fits
import numpy

from ixpeobssim.binning.base import xEventBinningBase, xBinnedFileBase
from ixpeobssim.binning.fmt import xBinTableHDULC, xBinTableHDUPP
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.core.hist import xScatterPlot
from ixpeobssim.utils.logging_ import logger, abort


### pylint: disable=invalid-name, attribute-defined-outside-init, too-few-public-methods
### pylint: disable=no-member, arguments-differ


class xEventBinningCMAP(xEventBinningBase):

    """Class for CMAP binning.
    """

    INTENT = 'count map in sky coordinates'
    SUPPORTED_KWARGS = ['mc'] + xEventBinningBase._image_wcs_kwargs()

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        self.process_image_ref_kwargs()

    def bin_(self):
        """Overloaded method.
        """
        ra, dec = self.event_file.sky_position_data(self.get('mc'))
        wcs_ = self._build_image_wcs(**self.kwargs)
        x, y, binning = self._pixelize_skycoords(ra, dec, wcs_)
        data, _, _ = numpy.histogram2d(x, y, bins=binning)
        header = wcs_.to_header()
        hdu_list = fits.HDUList([self.build_primary_hdu(data, header)])
        self.write_output_file(hdu_list)



class xBinnedMap(xBinnedFileBase):

    """Display interface to binned CMAP files.

    While this is essentially an xFITSImage, we need to inherit from
    xBinnedFileBase, as the latter provides the implementation of the
    from_file_list() slot. Instead of inheriting from xFITSImage, too, we
    preferred composition, instead.

    Warning
    -------
    I am more and more convinced that this is poor design, and maybe the pylint
    super-init-not-called is there for a reason. We should probably refactor the
    class and make the implementation cleaner.
    """

    # pylint: disable=super-init-not-called, abstract-method
    def __init__(self, file_path):
        """Constructor.
        """
        self.fits_image = xFITSImageBase(file_path)
        # These are needed to equip the class instance with all the members
        # that are necessary to call the wite() method of xBinnedFile.
        self.__data_dict = {}
        self._xBinnedFileBase__data_dict = self.__data_dict
        self.hdu_list = self.fits_image.hdu_list

    def __iadd__(self, other):
        """Overloaded method for CMAP binned data addition.
        """
        self._check_iadd(other)
        self.fits_image += other.fits_image
        return self

    def plot(self, **kwargs):
        """Plot the data. The kwargs passed to plt.imshow().
        """
        kwargs.setdefault('zlabel', 'Counts/pixel')
        return self.fits_image.plot(**kwargs)



class xEventBinningLC(xEventBinningBase):

    """Class for LC binning.
    """

    INTENT = 'light curve'
    SUPPORTED_KWARGS = ['tmin', 'tmax', 'tbins', 'tbinalg', 'tbinfile']

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        if self.get('tmin') is None:
            self.set('tmin', self.event_file.start_met())
        if self.get('tmax') is None:
            self.set('tmax', self.event_file.stop_met())

    def make_binning(self):
        """Build the light-curve binning.
        """
        kwargs = dict(min_val=self.get('tmin'), max_val=self.get('tmax'),
                      num_bins=self.get('tbins'), bin_file=self.get('tbinfile'))
        return xEventBinningBase.make_binning(self.get('tbinalg'), **kwargs)

    def bin_(self):
        """Overloaded method.
        """
        counts, edges = numpy.histogram(self.event_file.time_data(),
                                        bins=self.make_binning())
        primary_hdu = self.build_primary_hdu()
        data = [self.bin_centers(edges), self.bin_widths(edges), counts, numpy.sqrt(counts)]
        rate_hdu = xBinTableHDULC(data)
        rate_hdu.setup_header(self.event_file.primary_keywords())
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, rate_hdu, gti_hdu])
        self.write_output_file(hdu_list)



class xBinnedLightCurve(xBinnedFileBase):

    """Binned light curve.
    """

    def _read_data(self):
        """Overloaded method.
        """
        self._read_binary_table_data(xBinTableHDULC.NAME)

    def __iadd__(self, other):
        """ Overloaded method for LC binned file addition.
        """
        self._check_iadd(other, ('COUNTS', 'ERROR'), ('TIME', 'TIMEDEL'))
        self.COUNTS += other.COUNTS
        self.ERROR = numpy.sqrt(self.ERROR**2. + other.ERROR**2.)
        return self

    def plot(self):
        """Overloaded plot method.
        """
        rate = self.COUNTS / self.TIMEDEL
        rate_err = self.ERROR / self.TIMEDEL
        xScatterPlot(self.TIME, rate, rate_err, xlabel='MET [s]', ylabel='Rate [Hz]').plot()



class xEventBinningPP(xEventBinningBase):

    """Class for pulse-profile binning.
    """

    INTENT = 'pulse profile'
    SUPPORTED_KWARGS = ['phasebins']

    def make_binning(self):
        """Build the light-curve binning.
        """
        kwargs = dict(min_val=0., max_val=1., num_bins=self.get('phasebins'))
        return xEventBinningBase.make_binning('LIN', **kwargs)

    def bin_(self):
        """Overloaded method.
        """
        phase = self.event_file.phase_data()
        if phase is None:
            logger.error('PHASE column not found, do you need to run xpphase?')
            abort('Cannot create pulse profile')
        counts, edges = numpy.histogram(phase, bins=self.make_binning())
        primary_hdu = self.build_primary_hdu()
        data = [self.bin_centers(edges), self.bin_widths(edges), counts,
                numpy.sqrt(counts)]
        rate_hdu = xBinTableHDUPP(data)
        rate_hdu.setup_header(self.event_file.primary_keywords())
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, rate_hdu, gti_hdu])
        self.write_output_file(hdu_list)



class xBinnedPulseProfile(xBinnedFileBase):

    """Binned pulse profile.
    """

    def _read_data(self):
        """Overloaded method.
        """
        self._read_binary_table_data(xBinTableHDUPP.NAME)

    def __iadd__(self, other):
        """ Overloaded method for PP binned file addition.
        """
        self._check_iadd(other, ('COUNTS', 'ERROR'), ('PHASE', 'PHASEDEL'))
        self.COUNTS += other.COUNTS
        self.ERROR = numpy.sqrt(self.ERROR**2. + other.ERROR**2.)
        return self

    def plot(self, **kwargs):
        """Overloaded plot method.
        """
        xScatterPlot(self.PHASE, self.COUNTS, self.ERROR, xlabel='Pulse phase',
                     ylabel='Counts').plot(**kwargs)

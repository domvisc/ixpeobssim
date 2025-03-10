.. _release_notes:
Release notes
=============


*ixpeobssim (29.0.0) - Thu, 08 Sep 2022 11:19:50 +0200*

First version for public release.


*ixpeobssim (28.4.0) - Sat, 27 Aug 2022 10:59:55 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/271
* New linearized error propagation in the polarization cube subtraction.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/614


*ixpeobssim (28.3.0) - Fri, 26 Aug 2022 12:43:50 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/270
* Structure refactored to facilitate splitting the IRF generation code out of
  ixpeobssim.
* New irf.ebounds module added, with the energy grid definition.
* irfgen import removed from the argparse_ module.
* GPD filling temperature and pressure moved to irfgen.gpd
* Names and path for the MMA effective area and vignetting files moved to irfgen.mma
* UV filter file naming moved to irfgen.du
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/612


*ixpeobssim (28.2.0) - Wed, 24 Aug 2022 07:10:03 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/269
* Creation of weighted polarization cubes, maps and map cubes with acceptance
  correction disengaged until we have proper arf files with the `SIMPLE`
  weighting scheme.
* New weighting_scheme() hook, defaulting to None, added to the xEffectiveArea class.
* xStokesAnalysis constructor signature changed, in preparation of the addition
  of the energy flux to the polarization cubes.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/613


*ixpeobssim (28.1.0) - Tue, 23 Aug 2022 18:37:19 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/268
* New infrastructure in place for generating and using 2d (as in "non azimuthally
  symmetric") PSF.


*ixpeobssim (28.0.1) - Fri, 19 Aug 2022 07:58:59 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/267
* Docs entry page revamped.
* As-run target list updated.


*ixpeobssim (28.0.0) - Thu, 28 Jul 2022 08:20:14 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/266
* Massive refactoring of the IRF files, updated to the latest structure and
  naming conventions.
* Obsolete response files removed altogether.
* Code resolving the paths for the response files largely simplified.
* caldb folder ixpe/mma moved to ixpe/xrt/bcf to match the real CALDB.
* Unit tests updated.
* IRF docs completely revised.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/564


*ixpeobssim (27.0.0) - Thu, 21 Jul 2022 16:31:29 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/263
* ``rotate`` keyword argument removed from all the functions in srcmodel.polarization
  (from now on, when passing input maps in angle or PD/PA components, it is assumed
  that the angle is measured from the celestial North).
* xy mode for reading in polarization maps and aligning Stokes parameters
  removed altogether---we're still accepting PD/PA, but we should really encourage
  people to work in normalized Q/U space for extended sources.
* Origin of coordinates for measuring the position angle now correctly set to
  the celestial North, for both the visualization and the Stokes alignment.
* New data structures in srcmodel.polarization for radial and tangential
  polarization fields with arbitrary radial profiles.
* New toy_radial_disk and toy_tangential_disk examples, illustrating the new
  functionality.
* ``casa` example renamed to ``toy_casa``, and fully revamped.
* Clocking direction of the DUs fixed.
* Additional 90 degree rotations added in the photon-list generation and in
  xpsimfmt to fix the orientation of the polarization patterns in the e2e
  workflow.
* Unit tests added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/349
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/596
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/597


*ixpeobssim (26.6.1) - Thu, 21 Jul 2022 07:11:18 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/264
* Small additions for the final version of the Software X paper.


*ixpeobssim (26.6.0) - Wed, 13 Jul 2022 12:37:11 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/257
* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/259
* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/261
* Photon list mechanism implemented for the xChandraObservation model component class.
* Vignetting now correctly applied in the photon list workflow.
* ``DETX`` and ``DETY`` columns added in the xpsimfmt output files.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/593
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/606
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/607


*ixpeobssim (26.5.0) - Wed, 13 Jul 2022 11:28:26 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/261
* Script to display the as-run target list added.
* New configuration file for an unpolarized point source.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/600


*ixpeobssim (26.4.0) - Wed, 13 Jul 2022 11:22:07 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/260
* Brute force workaround for a regression introduced in numpy 1.22.0
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/608


*ixpeobssim (26.3.3) - Mon, 30 May 2022 12:07:45 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/256
* Small fix in the polarization cube subtraction and multiplication (credits: Lawrence P.)


*ixpeobssim (26.3.2) - Wed, 18 May 2022 15:40:05 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/255
* Minor fix to the zlabel for binned count maps.


*ixpeobssim (26.3.1) - Wed, 18 May 2022 13:59:39 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/254
* Bug fix in xpphase for ephemeris referred to times before the start of the
  observation.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/601


*ixpeobssim (26.3.0) - Wed, 18 May 2022 13:29:57 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/252
* New ixpeobssim.event.evt.xEventFileFriend class added to handle synchronized pairs
  of level-1 and level-2 data.
* xpselect generalized to accept a binary selection mask.


*ixpeobssim (26.2.0) - Wed, 18 May 2022 12:27:15 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/250
* Polarization map cubes equipped with a general-purpose convolution routine and
  a plotting hook for the significance.
* Circular sum kernel facility added.
* Alignment and radial profile of polarization maps for polarization map cubes.


*ixpeobssim (26.1.1) - Thu, 12 May 2022 17:28:16 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/253
* Fix for issue 599
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/599


*ixpeobssim (26.1.0) - Sat, 07 May 2022 21:07:09 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/249
* (This is an intermediate release with the specific purpose of keeping track
  of the analysis for the magnetar discovery paper, and we are not closing any
  of the related issues, just yet.)
* Added fiducial backscal value.
* Inhibit circle/annuli and ds9 region file selections, and writing the BACKSCAL
  header keyword for the first.
* Propagating the BACKSCAL value to the binned polarization cubes.
* Initial implementation of the PCUBE subtraction.
* Script for the 4u analysis for the Science paper.


*ixpeobssim (26.0.1) - Mon, 02 May 2022 16:06:51 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/248/
* Small test script added.


*ixpeobssim (26.0.0) - Mon, 02 May 2022 15:45:08 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/247/
* Default IRF name bumped to ixpe:obssim:v11
* PSF parametrization changed: maximum radius pushed out to 480 arcseconds,
  manual scale factors removed, and parametrizations set for the three MMAs
  separately.
* Docs updated.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/158
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/580
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/428


*ixpeobssim (25.8.0) - Thu, 28 Apr 2022 09:25:27 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/246
* Deterministic implementation of xppicorr.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/595


*ixpeobssim (25.7.0) - Fri, 22 Apr 2022 19:15:16 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/245
* Supporting regions version 0.6.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/589


*ixpeobssim (25.6.5) - Fri, 01 Apr 2022 12:14:21 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/244
* Specific target for testing the local installation added to the Makefile.
* Minor bug fix.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/590


*ixpeobssim (25.6.4) - Fri, 01 Apr 2022 12:10:16 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/242
* Rendering of the LTP recast in terms of the TWGs.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/578


*ixpeobssim (25.6.3) - Wed, 09 Mar 2022 15:36:06 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/243
* Removing all time-dependent PI correction files, as they now live in a separate
  repository: https://bitbucket.org/ixpesw/pi_corr_caldb/
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/579


*ixpeobssim (25.6.2) - Mon, 07 Mar 2022 15:09:52 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/241
* Method to calculate the energy spectrum added to the magnetar model interface.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/388


*ixpeobssim (25.6.1) - Tue, 01 Mar 2022 10:09:07 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/240
* Minor changes.


*ixpeobssim (25.6.0) - Mon, 28 Feb 2022 14:43:32 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/239
* New xpstokesrandom and xpstokesshuffle implemented, and added to the docs and
  wrapped in the pipeline.
* Unit tests added.


*ixpeobssim (25.5.0) - Mon, 28 Feb 2022 14:13:39 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/238
* New plotting style for the polarization cubes (the API should be considered
  experimental and might evolve as we learn to use the new features).
* New features in the ixpeobssim.utils.matplotlib_ module to support the new
  plotting of polarization cubes.


*ixpeobssim (25.4.0) - Mon, 28 Feb 2022 13:59:49 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/237
* Proper calculation of the detection significance added to polarization cubes,
  maps and map cubes.
* Additional fields P_VALUE and CONFID added to the polarization cubes, maps
  and map cubes.
* A few fix for zero-division errors.
* Small fix in summing the values of N_EFF and FRAC_W across polarization cubes,
  maps and map cubes.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/467


*ixpeobssim (25.3.4) - Fri, 25 Feb 2022 14:15:43 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/236
* New PI correction files for the ~complete Cas A observation 01001301


*ixpeobssim (25.3.3) - Fri, 25 Feb 2022 09:36:45 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/235
* Fix for issue #574.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/574


*ixpeobssim (25.3.2) - Thu, 24 Feb 2022 06:20:16 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/234
* (Yet another) fix for bug #567.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/567


*ixpeobssim (25.3.1) - Wed, 23 Feb 2022 14:42:40 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/233
* De-correction for the GPD quantum efficiency applying when simulating photon lists
  for the instrumental background---see, however all the caveats in the
  relevant issue.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/517


*ixpeobssim (25.3.0) - Wed, 23 Feb 2022 14:04:03 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/231
* Full refactoring of the code dealing with the response files.
* All format specifications moved to ixpeobssim.irfgen.fmt
* New xSpecRespBase class added, acting as a base class for the effective area,
  modulation factor and modulation response function, and equipped to use the
  SYS_MIN and SYS_MAX columns, when available.
* Vignetting factored out of the effective area class.
* pylinted.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/568
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/313


*ixpeobssim (25.2.1) - Tue, 22 Feb 2022 16:15:04 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/232
* Bug fix for issue #567
* Phase and time grids in the photon list now driven by the proper class members.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/567


*ixpeobssim (25.2.0) - Fri, 18 Feb 2022 11:35:38 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/230
* xppiscale.py renamed as xppicorr.py and generalized to global scale and offset
  corrections, as well as generic time-dependent corrections driven from a
  FITS file.
* Initial PI correction for the first chunck of Cas A observation (v02) added.
* Unit tests added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/568


*ixpeobssim (25.1.0) - Fri, 18 Feb 2022 10:36:58 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/229
* New facilities related to the exposure calculation.
* Livetime cube binning algorithm (LTCUBE) added---it saves a ltcube with information
  on the elapsed time in each theta bin in each pixel of the map.
* LTCUBE supported in xpbinview
* New xpancrkey and xpexposure apps added.
* New toy_offaxis configuration file and associated analysis pipeline illustrating
  the new exposure functionality.
* IN_SAA and TARGET_OCCULT columns in the SC_DATA extensions now driven by
  whether we are taking data---they are identically zero if the --saa and/or the
  --occult flags are set to False from command line.
* Bug fix in the __iadd__() slot for xBinnedAreaRateMap objects.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/433
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/572
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/562


*ixpeobssim (25.0.0) - Thu, 17 Feb 2022 13:46:16 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/223
* Adding support for weights in binned MDP and polarization maps.
* Added --acceptcorr option in xpbin for polarization data products.
* Small tweak to the binning functions to accept tuples in addition to lists.
* Fix for issue #443.
* Modulation cubes and associated classes removed.
* Major change to the format of the polarization cubes and the MDP and
  polarization maps and map cubes to keep track of all the necessary figures
  for a correct calculation of the MPD, as well as for holding errors on the
  Stokes parameters and the significance of a polarization measurement.
* Supporting errors on Stokes parameters in polarization cubes and maps.
* Bug fix in xpevtstat.py
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/443
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/540
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/565
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/566
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/444


*ixpeobssim (24.0.0) - Tue, 08 Feb 2022 15:15:19 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/227
* Massive cleanup of the ixpeobssim.evt.event module, with parts moved out to the
  new ixpeobssim.evt.fmt and ixpeobssim.evt.gti modules, and a few obsolete
  interfaces, such as _radec_to_xy_int(), removed.
* Major cleanup of the xEventList class, with obsolete interfaces removed, and
  WCS information used consistently throughout.
* All WCS-related header keywords are now consistently set through the proper
  keyword arguments of the fits.Column objects, as opposed to manually hacking
  the header itself.
* WCS information added to the output xpsimfmt files, that should be now
  properly displayed in ds9.
* Comprehensive revision of the binary table headers for simulated files.
* ixpeobssim.evt.ixpesim streamlined.
* build_wcs() signature changed for consistency.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/518
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/523
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/526
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/538
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/548
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/550
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/552


*ixpeobssim (23.8.1) - Tue, 08 Feb 2022 11:32:43 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/228
* Documentation on the binned data products fully revamped.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/478
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/551


*ixpeobssim (23.8.0) - Tue, 08 Feb 2022 10:15:04 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/226
* Bug fix for the phase folding returning values outside the interval [0, 1]
* Avoid applying the vignetting twice in the Chandra-to-IXPE converter.
* Modified Chandra-to-IXPE workflow using the exposure.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/131
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/488


*ixpeobssim (23.7.0) - Fri, 04 Feb 2022 16:09:45 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/225
* Added granular invert command line switches to xpselect; this allow to
  invert (i.e., take the logical not) of any of the selections applied.
* Small bug fix.
* Unit tests added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/549


*ixpeobssim (23.6.1) - Thu, 03 Feb 2022 10:26:54 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/224
* Minor change to toy_pollin to match the polling definition in XSPEC.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/547


*ixpeobssim (23.6.0) - Wed, 02 Feb 2022 17:10:33 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/221
* energy_spectrum changed to photon_spectrum throughout.
* Index for the EXB changed.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/544


*ixpeobssim (23.5.0) - Wed, 02 Feb 2022 16:56:24 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/222
* Small refactoring in the binning routines, with all the I/O dictionaries
  moved into the binning module, so that they can be effectively used in the apps.
* Mechanism for building the path to the output file fixed for consistency with
  the other apps in xpphase.py and xpophase.py
* New set_column() class member added to evt.event.xEventFile.
* New app xppiscale.py added, and included in the pipeline facilities.
* New toy_ms_pulsar configuration file and associated pipeline added.
* Small tweak to the header keywords in xpsimfmt.py.
* xpphotonlist added to the pipeline facilities.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/439


*ixpeobssim (23.4.0) - Wed, 02 Feb 2022 08:34:21 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/219
* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/220
* Bug fix for filtering mismatch in weighted polarization analysis.
* Minor tweaks to the rendering of bivariate splines.
* Added a small macro with preliminary plots for the ixpeobssim paper.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/541


*ixpeobssim (23.3.0) - Fri, 28 Jan 2022 15:03:30 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/218
* xStokesAnalysis modified to filter out malformed events and events outside
  the 0--15 keV energy where the response functions can be sensibly extrapolated.
* Improved diagnostics in xStokesAnalysis.
* More sensible error message from xpbin.py when the input file is not found.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/539


*ixpeobssim (23.2.1) - Wed, 26 Jan 2022 18:49:36 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/217
* Improved errore messages when failing consistency check in summing binned products.
* Using TSTART and TSTOP as default values for LC binning bounds (as opposed to
  the extremes of the GTIs, which are generally different for the three detectors
  in the same observation).
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/537


*ixpeobssim (23.2.0) - Wed, 26 Jan 2022 17:04:18 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/216
* New background PHA1 files created from the first observation of SMC X-1, and
  code added to perform the proper scaling to create usable templates.
* New xTemplateInstrumentalBkg class for generating template-driven background spectra.
* Added (and enabled by default) an option to prevent the convolution with the
  instrumental background spectrum with the energy dispersion, and modified the
  handling of the energy bounds for the simulation in the two cases.
* Docs added.
* Sample configuration file (instrumental_bkg_smcx1) added to illustrate the new
  functionality.
* Realistic instrumental background added to the Cas A configuration file.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/535


*ixpeobssim (23.1.0) - Wed, 26 Jan 2022 10:30:25 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/215
* ixpeobssim.evt.subselect refactored to interoperate with filtered, level-2 data.
* Livetime correction disabled by default in xpselect.
* Docs and unit tests updated.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/536


*ixpeobssim (23.0.1) - Tue, 25 Jan 2022 18:49:28 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/214
* LTP updated.
* Last update label added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/533


*ixpeobssim (23.0.0) - Mon, 24 Jan 2022 13:34:04 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/213
* This is the first release capable of operating on flight data, and most of the
  changes originate from the very first experience with the Cas A data.
* Using OBJ_RA and OBJ_DEC (rather than PNT_RA and PNT_DEC, that are not present
  in the level-2 file primary header) as a default value for centering WCS objects.
* Kislat analysis re-cast in Q and U, as opposed to PHI.
* EXPOSURE keyword removed from the event lists, and added at the xpbin.py level
  to allow for fitting in XSPEC.
* Pixel grid definition for X and Y changed from 900 x 900 pixels at 2 arcsec steps
  to 600 x 600 pixels at 2.6 arcsec steps.
* Physical energy in keV retrieved via the PI column in event files.
* RA and DEC coordinates retrieved via X and Y in event files.
* xpselect refactoring to handle with the fact that the LIVETIME columns is not
  included in filtered level-2 event lists.
* Minor changes.
* Docs updated.
* Unit tests added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/531
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/528
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/532
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/530
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/529


*ixpeobssim (22.0.0) - Sun, 23 Jan 2022 09:21:26 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/210
* Full refactoring of ixpeobssim.evt.subselect in order to correctly propagate
  the livetime through the time and phase selections.
* xEventFile.average_deadtime_per_event() hook added.
* --phimin, --phimax and --invert options removed from xpselect.
* Added livetime-correction options to xpselect.
* Headers properly updated in xpselect.
* Time-related keywords added to the MONTE_CARLO extension.
* Infrastructure to control the count spectrum spline (ny, kx and ky) in place.
* New xStepFunction class added.
* New livetime examples (with selection in time and phase) revised.
* Documentation section about xpselect added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/378
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/55
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/514
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/169


*ixpeobssim (21.4.0) - Sun, 23 Jan 2022 07:46:55 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/212
* New xpstokessmear application added to test the effect of the spurious modulation
  correction via a gaussian smearing of the Stokes parameters.
* New xpaddmofweights application added to process a level-2 file adding a new
  column with weights based on the modulation factor as a function of the energy.
* Output support enhanced in the xEventList class, via the addition of the
  add_column(), add_columns(), remove_columns() and write() methods.
* General binary search method to locate bin and bin values in multi-dimensional
  histograms added.
* xpstokessmear and xpaddmodfweights added to the reference docs.
* xpaddmofweights and xpstokessmear added to the pipeline.
* Minor refactoring of the basic app structure.
* Unit tests added.
* Copyright notice updated.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/512
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/521


*ixpeobssim (21.3.3) - Thu, 20 Jan 2022 22:02:12 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/211
* Minor doc update.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/506


*ixpeobssim (21.3.2) - Wed, 19 Jan 2022 20:39:07 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/209
* Bug fix for having xpsimfmt inter-operate with event files with no MONTE_CARLO
  extension.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/515


*ixpeobssim (21.3.1) - Tue, 18 Jan 2022 05:53:58 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/208
* Emergency patch for issue #513 (photon-list mechanism not working with instrumental
  background).
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/513


*ixpeobssim (21.3.0) - Mon, 17 Jan 2022 18:39:51 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/207
* Sorting the photon list before writing them out to FITS, fixing a fairly
  serious flaw in the mechanism.
* Roll angle added to the SC_DATA extension.
* SC_DATA extension added to the photon lists.
* Adding RA, DEC, X and Y to the xpsimfmt output file.
* Properly handling dithering, vignetting and fiducial cut in the photon lists.
* Using scdata=False for the test_instrumental_background test.
* pointing_ra/dec changed to ra/dec_pnt throughout.
* Polarization angle in the photon lists rotated to the GPD reference frame, and
  inverse transformation implemented in xpsimfmt
* DET_ID overwritten by xpsimfmt
* Added option to use MC/reconstructed absorption point in xpsimfmt.
* Documentation updated and unit test added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/494
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/498
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/500



*ixpeobssim (21.2.0) - Sat, 15 Jan 2022 10:48:30 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/203
* Implemented dithering directly to the pointing direction, so that it gets
  propagated to the SC_DATA binary table.
* Vignetting now correctly applied.
* Moved parse_dithering_kwargs() method to the ixpeobssim.instrument.mma module.
* Command-line options refactored.
* Added facility to recover the pointing direction from the SC_DATA table.
* LAUNCH_DATE and LAUNCH_MET added in the time_ module.
* Added unit tests.
* Added documentation.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/431
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/432


*ixpeobssim (21.1.2) - Thu, 13 Jan 2022 16:34:48 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/206
* Bugged commit https://bitbucket.org/ixpesw/ixpeobssim/commits/264dad9b5b1549ec83d9a2dfb874491ee3901045
  reverted.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/487
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/497


*ixpeobssim (21.1.1) - Thu, 13 Jan 2022 15:19:10 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/205
* XSPEC version parsing removed.
* Added unit test.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/507


*ixpeobssim (21.1.0) - Thu, 13 Jan 2022 14:59:08 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/204
* Default suffix for xpstokesalign changed from '_phialign' to '_stokesalign'
* Preventing xpstokesalign from changing the DEPHI column, if present in the input
  event list.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/504
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/505


*ixpeobssim (21.0.1) - Tue, 11 Jan 2022 20:55:54 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/201
* Emergency fix for the setuptools files after the refactoring of the CALDB.
* ixpeobssim/srcmodel/par_files moved to ixpeobssim/srcmodel/parfiles.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/499


*ixpeobssim (21.0.0) - Mon, 10 Jan 2022 17:49:05 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/201
* Full re-organization of the pseudo CALDB to match the structure of the real CALDB.
* IRF-name separator changed from "_" to ":" to allow CALDB-like file names and
  properly support weights.
* File name conventions for the IRF files aligned with the real CALDB starting
  from version 10.
* New keywords (and checksum) added to all the response files.
* COMMENT fields pertaining to the version and weight for response functions
  removed from all the headers, since this information is now tracked in proper
  keywords.
* IRF documentation fully revamped.
* Default IRF name set to "ixpe:obssim:v10".
* Modulation response function added to the xIrfSet class.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/496
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/462
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/468
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/479
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/492


*ixpeobssim (20.2.0) - Fri, 07 Jan 2022 16:23:47 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/200
* Added facilities to display the observation plan.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/495


*ixpeobssim (20.1.0) - Fri, 07 Jan 2022 16:16:00 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/199
* Polarization alignment according to an input model re-casted in terms of the
  Stokes parameters.
* xpphialign.py renamed as xpstokesalign.py
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/493


*ixpeobssim (20.0.0) - Fri, 17 Dec 2021 16:02:54 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/198
* All multiplicative polarization models renamed and in synch with the XSPEC repository:
      * ``constpol`` is now ``polconst``
      * ``linpol`` is now ``pollin``
      * ``powpol`` is now ``polpow``
      * ``quadpol`` model has bee removed.
      * (Note that the parameter names have been changed, as well)
* `Do not forget to cleanup and recompile the ixpeobssim local models!`
* All model names changed in the codebase.
* Documentation updated.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/489


*ixpeobssim (19.7.0) - Fri, 17 Dec 2021 10:18:46 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/197
* IXPE TLE updated with the first post-launch values.
* Docs updated.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/490


*ixpeobssim (19.6.1) - Fri, 17 Dec 2021 08:49:42 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/196
* Minor fixes to the documentation.


*ixpeobssim (19.6.0) - Tue, 30 Nov 2021 10:10:23 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/195
* Initial implementation of the animation module.
* xpsonify improved to support animation.
* Docs updated.


*ixpeobssim (19.5.3) - Tue, 30 Nov 2021 10:06:23 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/194
* Fix a runtime zero-division error in xpbin.py
* Fix a runtime error in core.fitsio.py
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/470
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/483


*ixpeobssim (19.5.2) - Thu, 25 Nov 2021 18:17:54 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/193
* Docs tweaked.


*ixpeobssim (19.5.1) - Wed, 24 Nov 2021 15:58:08 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/192
* Bug fix.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/482


*ixpeobssim (19.5.0) - Wed, 24 Nov 2021 13:44:05 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/191
* New sonification module.
* New application to transform a photon list into an MIDI and/or audio file.
* Documentation updated.


*ixpeobssim (19.4.1) - Mon, 22 Nov 2021 18:39:23 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/190
* Unit test added for issue #179 (invalid)
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/179


*ixpeobssim (19.4.0) - Mon, 22 Nov 2021 16:26:09 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/189
* xInstrumentalBackground class refactored to support photon lists.
* Energy bounds for the instrumental background now correctly inferred from the
  simulation setup.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/481


*ixpeobssim (19.3.3) - Mon, 22 Nov 2021 11:50:31 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/188
* Toy response functions removed.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/452


*ixpeobssim (19.3.2) - Mon, 22 Nov 2021 11:29:01 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/187
* main() entry point added to xpphotonlist
* Unit test added to ensure that all the apps have appropriate entry points to
  run in user mode.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/480


*ixpeobssim (19.3.1) - Mon, 22 Nov 2021 09:29:03 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/186
* A couple of typos fixed.
* Energy spectrum changed to photon spectrum throughout.
* Docs on binary systems added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/430
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/261
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/384
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/386


*ixpeobssim (19.3.0) - Sun, 21 Nov 2021 20:37:10 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/185
* highecut_power_law spectral model added.
* xpchrgmap added to the pipeline, and corresponding command-line parser modified accordingly.
* Simplified observing plan simulation implemented.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/449
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/450


*ixpeobssim (19.2.5) - Sat, 20 Nov 2021 17:21:56 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/184
* XSPEC headers moved to a separate file to facilitate supporting multiple XSPEC version.
* Collecting PyXspec and XSPEC version strings.
* Conditional compilation for the XSPEC headers, to support the new include layout
  in XSPEC version 12.12.0.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/472


*ixpeobssim (19.2.2) - Sat, 20 Nov 2021 09:05:34 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/183
* Fix for malformed TLE in sgp4 version 2.20
* TLE epoch changed from January 1, 2021 to December 9, 2021.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/473


*ixpeobssim (19.2.1) - Fri, 19 Nov 2021 11:19:28 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/182
* Minor fixes to the docs.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/476


*ixpeobssim (19.2.0) - Fri, 19 Nov 2021 09:50:30 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/180
* xpchrgcorr.py removed in favor of the official tool available in gpdsw.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/475


*ixpeobssim (19.1.0) - Fri, 19 Nov 2021 09:32:43 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/181
* CI Python version changed from 3.9.6 to 3.6.15.
* A couple of tweaks to support Python 3.6.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/477


*ixpeobssim (19.0.0) - Thu, 18 Nov 2021 13:54:55 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/179
* First implementation of the photon list mechanism.
* xBaseEventList class added, and xEventList refactored to support the new xPhotonList.
* Keeping track of the primary header comments in the IRF files.
* xEventList.filled_array() method removed.
* Refactoring of the roi module.
* Docs updated.


*ixpeobssim (18.1.1) - Thu, 18 Nov 2021 10:56:19 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/178
* main() entry point added for consistency to all the apps, see
  https://bitbucket.org/ixpesw/ixpeobssim/issues/469


*ixpeobssim (18.1.0) - Wed, 17 Nov 2021 14:14:14 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/176
* Bug fix in the flux normalization for the magnetar models.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/471


*ixpeobssim (18.0.1) - Wed, 17 Nov 2021 13:22:41 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/176
* Bug fix in the xpbin.py pixsize command-line switch.
* --dpi option added to xpbinview.py

* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/454


*ixpeobssim (18.0.0) - Wed, 13 Oct 2021 15:59:48 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/175
* W_MOM column added to the event lists (provisionally set to 1).
* DET_Q and DET_U changed into Q and U.
* X and Y changed from int to floats.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/415
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/424


*ixpeobssim (17.3.0) - Fri, 08 Oct 2021 09:29:13 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/174
* New MMA effective area curves with a refined analysis of the MMA calibration data.
* IRF v9 created.
* Version number added to the IRF heders as a comment.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/460


*ixpeobssim (17.2.0) - Thu, 07 Oct 2021 14:49:07 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/173
* xpbin.py generalized to support weights.
* xpsimfmt.py generalized to support weights and added to the pipeline, with some
  tweaks to allow for a full spectro-polarimetric fit in XSPEC.
* xpcustomirf.py generalized to support weights, and now generating a nominal
  vignetting function to allow the loading of the effective area.
* Bug fix in xpevtstat.py when running on input files with no MONTE_CARLO extension.
* v8 iteration of the response functions added (but not the default, yet). Note
  this is the version passed over to the SOC to start populating the CALDB, and
  includes the first set of response functions with weights.
* DET_ID added to the primary header of the IRF files.
* Keeping track of N_EFF and FRAC_W in the Stokes analysis a la Kislat.
* Major restructuring of the irfgen code to support the generation of response
  functions with weights.
* Small fix for the modified XSPEC errors.
* N_EFF and FRAC_W columns added in the polarization and modulation cubes.
* All references to the standard cut efficiency removed.
* Secular pressure values updated.
* Minor tweak to the ``utils.argparse_`` module.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/464
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/463
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/459


*ixpeobssim (17.1.0) - Thu, 02 Sep 2021 12:49:35 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/172
* New folder ``obsdata`` added to the hierarchy to hold observation-specific
  files (e.g., charging maps)
* Added vanilla charging maps, with all the values set to zero, to be used
  in observations where the detectors are initially completely discharged
  (and to be used by default).
* Charging parameters now read from the proper file in the preudo-CALDB.
* Charging-specific command-line switches modified (note that chrgtsteps is
  now chrgtstep, and we're setting the width of the step, rather than the
  number of steps).


*ixpeobssim (17.0.0) - Wed, 01 Sep 2021 15:06:39 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/171
* GTI calculation completely refactored.
* OCTI calculation revamped, now inhibiting the activation of the calibration
  sources in the SAA.
* New ``TIMELINE``, ``SC_DATA`` and ``OCTI`` extensions added to the output files
  to keep track of the status of the instrument along the orbit.
* New xpobsview application added for a quick look of a given observation timeline.
* xpobssim command-line switched tweaked for consistency, and new switches to
  control the minimum duration and the padding of the GTIs ans the OCTIs added.
* Livetime-related keywords fixed when the on-orbit calibration sources are
  activated (issue 457).
* Docs updated.
* Data format documentation regenerated as part of the docs creation (issue #429).
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/429
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/409
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/417
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/425
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/457


*ixpeobssim (16.17.0) - Tue, 24 Aug 2021 07:56:22 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/170
* Confidence interval calculation added in XSPEC fitting, and enabled by default.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/346


*ixpeobssim (16.16.0) - Thu, 19 Aug 2021 10:09:20 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/169
* Models for all the data challenge 1 sources added.
* Source documentation added.
* New gaussian line spectral model.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/410


*ixpeobssim (16.15.0) - Sat, 14 Aug 2021 19:14:23 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/168
* Added an optional ‘side’ argument to the bisect function in hist.py, matching
  the signature of numpy.searchsorted (default is ‘left’, as in numpy.searchsorted,
  so the change is backward-compatible).
* Small change in how events are assigned to the correct gain value by the gain()
  function in charging.py, to match the fact that the the self.__gain_data attribute
  now has the dimension of its time axis increased by one, matching exactly the
  time binning of the energy flux cube.
* Implemented the slow charging process---for now its parameters are hard-coded
  to zero, so that only the fast charging is actually active. We will fully enable
  the slow charging process when charging parameters will be taken from a CALDB
  file, see
  https://bitbucket.org/ixpesw/ixpeobssim/issues/449/set-the-parameters-for-the-charging-model
* Writing the slow charging map to the CHRG_MAP extension.


*ixpeobssim (16.14.0) - Sat, 14 Aug 2021 11:32:29 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/164
* Added an option to provide a list of input charging maps to xpobssim.py,
  along with the facilities to parse them.
* Modified the charging model in charging.py to accept the input charging map
* Added in charging.py two classes representing the PRIMARY and CHRG_MAP extension
  of the FITS charging map files
* Modified most of the functions in ``utils.time_`` to optionally accept a custom format
  (defaulting to DATETIME_FMT).
* Small app added for extracting charging maps from observation files and save
  them in a dedicated file.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/445


*ixpeobssim (16.13.0) - Fri, 13 Aug 2021 21:00:50 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/167
* Major restructuring of the auxiliary infrastructure for the response functions,
  but no change in any of the standard applications.
* Enhanced support for ixpeobssim-ixpesim inter-operation.
* AUX_VERSION bumped to version 3.
* New xpsimfmt.py and xpcustomirfs.py applications added.
* PI calculation for allx data sets improved, and bookkeeping added.
* Window contaminants correction implemented in xpsimspec.py
* Docs updated.


*ixpeobssim (16.12.1) - Mon, 09 Aug 2021 15:42:35 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/166
* A bunch of facilities related to spurious modulation added, with no
  change in any of the standard applications.


*ixpeobssim (16.12.0) - Wed, 04 Aug 2021 19:21:09 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/165
* Bug fix in handling magnetar tabular models.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/453


*ixpeobssim (16.11.0) - Thu, 15 Jul 2021 14:10:25 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/163
* xpphialign.py now changing phi, q and u consistently.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/441


*ixpeobssim (16.10.2) - Wed, 23 Jun 2021 18:30:23 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/159
* Draggable colorbar class added.
* Option for non-linear colorscale added in xpbinview.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/427


*ixpeobssim (16.10.1) - Wed, 23 Jun 2021 18:12:51 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/161
* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/162
* Fix for the generation of the magnetar model tables.
* Added docs for the argparse odd corner with negative number in engineering format.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/438


*ixpeobssim (16.10.0) - Thu, 10 Jun 2021 11:51:51 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/157
* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/158
* RA_PNT and DEC_PNT keywords used for the default ROI center in xpbin,
  xpselect and xpphialign
* aux.py renamed to auxiliary.py to allow interoperability with Windows.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/426
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/434


*ixpeobssim (16.9.1) - Wed, 09 Jun 2021 17:53:13 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/160
* Bug fix in how the model files were handle by the pipeline xpphialign wrapper.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/353


*ixpeobssim (16.9.0) - Thu, 03 Jun 2021 18:26:51 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/154
* Critical bug fix affecting the vignetting application, xpselect.py and the
  Chandra to IXPE conversion (please update).
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/423


*ixpeobssim (16.8.5) - Thu, 03 Jun 2021 17:50:50 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/155
* notebooks folder removed.
* Added a paragraph about the regions installation on the docs.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/418
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/405


*ixpeobssim (16.8.4) - Mon, 31 May 2021 08:19:33 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/153
* Normalization factor for the Galactic X-ray background changed.
* Interface to the ROSAT PSPC response matrix added.
* Unit test added.


*ixpeobssim (16.8.3) - Sat, 29 May 2021 10:17:24 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/152
* A series o minor tweaks to the data format, and in particular the header keywords.
* OBJECT, RA_PNT/RA_OBJ and DEC_PNT/DEC_OBJ keywords added, and xpobssim.py
  equipped with a new --objname command-line switch.
* DET_ID added for the physical identification of the detector units.
* Timing keywords updated in the GTI extension.
* APID, PKTTYPE and PKTSTYPE keyords removed.
* DAQ_VER keyword removed.
* CREAT_ID keyword removed, and version written into CREATOR.
* RUN_ID and STA_ID keywords removed.
* Header keyword comments capitalized.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/419
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/422


*ixpeobssim (16.8.2) - Sat, 29 May 2021 09:37:58 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/151
* New tool xpstat.py added for a quick look at the counts for a various components in
  a photon list.
* Docs updated.
* Minor refactoring of the energy binning base routine.
* Figure generation inhibited in a few unit tests.


*ixpeobssim (16.8.1) - Fri, 28 May 2021 14:16:29 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/150
* Minor tweaks to the IRF plotting.


*ixpeobssim (16.8.0) - Fri, 28 May 2021 14:12:28 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/149
* Initial implementation of the classes for the Extra-Galactic and Galactic
  X-ray background.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/57


*ixpeobssim (16.7.1) - Tue, 25 May 2021 20:17:28 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/148
* xpsimspec.py application added to create user spectrum files to be fed into ixpesim.
* Unit test for the energy redistribution added.


*ixpeobssim (16.7.0) - Mon, 24 May 2021 20:32:51 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/147
* New version of the IRF (v7) generated (but not the dafault, yet) with a
  non-diagonal response matrix, no 80% cut and Monte Carlo based modulation factor.
* Full machinery for processing and post-processing auxiliary files informing
  the response functions.
* Script to generate response functions at an arbitrary pressure added.
* New xLogNormal, xGeneralizedGaussian and xHat models added to core.modeling
* New xInterpolatedPiecewiseUnivariateSpline class added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/402
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/154


*ixpeobssim (16.6.1) - Tue, 18 May 2021 12:42:22 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/146
* Emergency patch for a regression in xpselect---this was a *MAJOR* breakage,
  if you have checked out version 16.0.0 please UPDATE IMMEDIATELY!
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/421


*ixpeobssim (16.6.0) - Thu, 13 May 2021 15:04:26 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/145
* New tool xpstripmc.py to process ixpeobssim photon lists and creating a
  verbatim copy without the Monte Carlo information (i.e., the MONTE_CARLO
  and ROI_TABLE extensions).
* --irfname option added to xpbin.py to support the analysis of files with no
  Monte Carlo information.
* xEventFile class modified to support photon lists with no Monte Carlo information.
* xpselect.py modified to support photon lists with no Monte Carlo information.
* Docs updated and unit tests added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/398
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/148


*ixpeobssim (16.5.0) - Tue, 13 Apr 2021 16:43:27 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/139
* Major refactoring of non-celestial sources, now split out the ixpeobssim.srcmodel.roi
  model into ixpeobssim.srcmodel.calibsrc
* First implementation of the FCW CalC source.
* Finalization of the event list refactored in its own method, automatically
  called right before the event list is written to file, to provide a unique
  and consistent interface for filling the ancillary columns.
* xpobssim modified to allow for interleaving celestial observations with FCW
  CalC calibration runs.
* New xpcalib.py app added to simulate calibration runs.
* Docs updates, unit tests added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/393
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/394
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/399


*ixpeobssim (16.4.1) - Tue, 13 Apr 2021 12:26:32 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/144
* Proper auxfile setup for the g21 and vela examples.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/408


*ixpeobssim (16.4.0) - Tue, 13 Apr 2021 11:18:46 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/141
* Dropping pyregion altogether in favor of the astropy affiliated package regions.
* Added sky filtering with astropy regions for sky coordinates.
* Docs, requirements and unit tests updated.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/381


*ixpeobssim (16.3.0) - Fri, 09 Apr 2021 16:06:33 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/142
* Changing PHE_Q and PHE_U columns to DET_Q and DET_U.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/400


*ixpeobssim (16.2.0) - Fri, 09 Apr 2021 16:03:39 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/143
* Configuration region files moved from ixpeobssim/config/fits to ixpeobssim/config/reg
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/233
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/406


*ixpeobssim (16.1.1) - Thu, 08 Apr 2021 08:29:04 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/140
* Fixed offset-by-one bug in the charging calculation.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/404


*ixpeobssim (16.1.0) - Tue, 06 Apr 2021 13:14:03 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/138
* Major refactoring of the code for generating response function, with lots of
  cleanup and complete removal of the toy response functions.
* Initial support for creating response functions from the full Monte Carlo
  simulation, with auxiliary files for the passive conversions and the
  energy dispersion, as well as the ROI size distribution.
* Major refactoring of the ixpeobssim.core.hist module, with a complete cleanup
  of the interfaces, support for errors on unweighted and weighted histogram and
  for data persistence in FITS format.
* KDE smoothing added to the histogram classes.


*ixpeobssim (16.0.0) - Thu, 01 Apr 2021 15:32:20 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/136
* New iteration (v6) of the response function---first one informed by the
  MMA and e2e calibration, and last one using the 80% cut.
* Use the MMA effective area curves from the mirror calibrations.
* Use the post-BAC best estimates of the GPD asymptotic pressures.
* Adjust the focal length to the measured value.
* Change the binning for the response functions.
* Use the measured PSF HPD
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/333
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/334
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/335
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/336
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/369
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/387


*ixpeobssim (15.1.0) - Wed, 31 Mar 2021 14:09:48 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/137
* Fix for the polarization degree negative values from the magnetar table models.
* Added support for magnetar models with QED off.


*ixpeobssim (15.0.0) - Tue, 23 Mar 2021 19:42:56 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/135
* NUM_CLU and LIVETIME columns added to the EVENTS extension.
* FILE_LVL keayword added to the file headers.
* Old livetime correction based on the number of discarded events replaced with
  the sum of event livetimes.
* A few methods related to the livetime added to the event file read interface.
* Small refactoring to avoid multiple conversions from start_date to start_met.
* pyregion import protected.
* Support for pseudo-Lv1a output added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/392


*ixpeobssim (14.2.1) - Tue, 23 Mar 2021 13:48:50 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/134
* Got rid of a few deprecation warnings from matplotlib 3.3
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/373


*ixpeobssim (14.2.0) - Mon, 22 Mar 2021 14:16:47 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/133
* Bug fix in the charging model with empty temporal bins.
* Bug fix in the charging model with a missing transpose.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/389


*ixpeobssim (14.1.0) - Mon, 22 Mar 2021 14:12:16 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/132
* Major refactoring of the xEphemeris class
* get_phase_func() removed
* phase_function() deprecated
* t0 argument to xpphase changed to met0
* xEphemeris.rvs() implemented, and unit test added.
* Ephemeris handling fixed in srcmodel.roi
* Periodic source examples cleaned up.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/52


*ixpeobssim (14.0.0) - Tue, 16 Mar 2021 11:26:05 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/131
* Added facilities to parse and operate with the magnetar models provided by
  Roberto and Roberto, and obsolete parsing routines removed.
* Added machinery for ixpeobssim auxiliary files.
* Example axp_1rxs_j1708.py revamped using the new functionality.
* Docs and unit tests updated.
* Some unintended fallout from merging pull request 129 cleaned up, and
  higher terms in the sourcse ephemeris disengaged until issue #52 is fixed.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/361


*ixpeobssim (13.1.0) - Fri, 12 Mar 2021 15:15:13 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/130
* New facility added for setting XSPEC model strings.


*ixpeobssim (13.0.0) - Sat, 27 Feb 2021 10:18:53 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/129
* Added the xpophase tool and xptimetophase replaced with xpphase.


*ixpeobssim (12.11.0) - Mon, 15 Feb 2021 17:04:03 +0100*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/128
* New xBinarySource class, and associated infrastructure and unit tests.
* New xpphase.py (replacing the old xptimetophase.py, now removed) and
  xpophase.py applications.
* New configuration file and associated pipeline toy_binary.py


*ixpeobssim (12.10.0) - Tue, 13 Oct 2020 13:31:28 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/127
* Added facility to save polarization map arrows as ds9 region file.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/361


*ixpeobssim (12.9.0) - Tue, 13 Oct 2020 09:20:09 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/125
* Added g21_bucciantini example and analysis pipeline.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/359


*ixpeobssim (12.8.0) - Mon, 12 Oct 2020 17:41:23 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/125
* New srcmodel.tdelays module and associated unit tests.


*ixpeobssim (12.7.0) - Sun, 04 Oct 2020 17:53:55 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/124
* Bug fix in xFITSImageBase.center(), now using the wcs conversions.
* New wcs facilities in utils.astro, and used in binning.
* Facilities to build intensity maps for arbitrary models in srcmodel.roi.
* Unit tests added.


*ixpeobssim (12.6.0) - Fri, 02 Oct 2020 15:54:42 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/98
* xEphemeris class moved out of srcmodel.roi into the new module srcmodel.ephemeris
* A bunch of facilities for binary sources added to srcmodel.ephemeris
* inverse() method implemented in the univariate base class.
* mjd_to_met() function added in ``utils.time_``
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/358


*ixpeobssim (12.5.1) - Thu, 01 Oct 2020 15:41:21 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/123
* New implementation of the angular separation function.


*ixpeobssim (12.5.0) - Thu, 01 Oct 2020 08:31:33 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/122
* Many improvements in xFITSImageBase plotting routines, courtesy of Niccolo Bucciantini.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/355


*ixpeobssim (12.4.0) - Thu, 24 Sep 2020 06:57:27 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/121
* New toy_pwn and toy_rim source examples added to aid the development of tools
  for the study of extended sources.
* Small tweaks to the binning module, and more zero-division-error protections added.
* evt.deconvolution module added (unused).
* xUniformAnnulus class added in srcmodel.roi
* Unit tests added.


*ixpeobssim (12.3.1) - Tue, 22 Sep 2020 10:24:26 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/120
* Bug fix in xFITSImageBase.sky_bounding_box(), courtesy of Niccolo Bucciantini.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/352


*ixpeobssim (12.3.0) - Sat, 19 Sep 2020 14:20:24 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/119
* Binned products can now be manipulated and saved to file.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/345


*ixpeobssim (12.2.0) - Sat, 19 Sep 2020 13:48:51 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/118
* Added the option to pass the tick marks labels on the colorbar of the
  xFITSImageBase class.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/351


*ixpeobssim (12.1.0) - Wed, 16 Sep 2020 06:40:30 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/117
* Smoothing out some rough edges around the XSPEC local models.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/350


*ixpeobssim (12.0.0) - Wed, 09 Sep 2020 15:02:09 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/116
* This is backward-incompatible release that incorporates significant
  changes and refactoring in several different areas.
* New IRF (version 5) released---incorporating a small change in the format
  of the modulation factor.
* Formalism in Kislat et al. (2015) now consistently implemented throughout.
* Major refactoring of the binned data structures: SCUBE algorithm removed,
  MCUBE algorithm deprecated, and several algorithms added (PHA1QN, PHA1N, PCUBE,
  MDPMAP, MDPMACUBE, PMAP, PMAPCUBE),
* Several small fixes in the binning routines.
* Improved support for XSPEC, including new models and purely polarimetric fits
  to the normalized Q/I and U/I Stokes parameters.
* Improved support for analysis and visualization of extended sources, including
  maps of MDP and normalized Stokes parameters.
* Small wrapper around GRPPHA added to the pipeline.
* Documentation updated.
* Documentation pdf target fixed.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/171
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/265
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/303
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/311
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/328
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/329
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/330
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/332
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/332
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/337
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/338
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/339
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/340
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/341
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/342
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/344


*ixpeobssim (11.2.1) - Thu, 20 Aug 2020 15:48:08 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/115
* Tentative fix for issue #325.
* Ephemeris for toy_periodic_source changed to trigger possible folding problems.
* Unit test added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/325


*ixpeobssim (11.2.0) - Thu, 20 Aug 2020 13:32:09 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/113
* New benchmark infrastructure added, with two examples.
* New xBinnedCountSpectrumSet implemented to calculate the broadband polarization,
  and releated facilities added.
* Docs updated.
* Command-line switch to initialize the fit parameters added back to xpxspec.py
* Added protection against wrong number of input files to xpxspec.py
* "pha1*" pattern now supported in pipeline.file_list(), and examples modified.
* Unit test streamlined.
* Resolution removed from the Gaussian model stat box.


*ixpeobssim (11.1.0) - Thu, 20 Aug 2020 07:39:34 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/114
* Write and read interfaces to MDP maps implemented.
* MDPMAP algorithm added to xpbin.py
* xpbinview.py now handling MPDPMAP binned files.
* Unit test added, and toy_disk example complemented.
* Docs updated.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/187


*ixpeobssim (11.0.0) - Tue, 18 Aug 2020 12:05:40 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/112
* Major rework of the XSPEC local models shipped with ixpeobssim, with the linpol
  and quadpol additions, and all parameter names changed for consistency.
* toy_linpol.py example added.
* Added facility to load the XSPEC local models programmatically.
* xpxspec and xpxspec3 merged and largely streamlined, with all the examples
  modified accordingly.
* pha1* tweak added to pipeline_file_list(), docs updated and examples modified.
* Docs for XSPEC support largely revised.


*ixpeobssim (10.5.1) - Wed, 12 Aug 2020 15:10:02 +0200*

* Ops---release notes updated.


*ixpeobssim (10.5.0) - Wed, 12 Aug 2020 15:08:12 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/111
* Cen A example cleaned up.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/310


*ixpeobssim (10.4.0) - Tue, 11 Aug 2020 15:12:47 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/109
* Specific versions added to requirements.txt
* Documentation updated with more details about the dependencies.
* Unit test added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/290


*ixpeobssim (10.3.0) - Tue, 11 Aug 2020 12:12:26 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/110
* Dropping the imp module in Python 3.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/277


*ixpeobssim (10.2.0) - Tue, 11 Aug 2020 12:05:50 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/107
* Support for scaling added in the xMDPRecord and xMDPTable classes.
* Broadband values for MDP tables now calculated dinamically at run time (this
  makes the bookkeeping much easier).
* eef and deadtime corrections added to xpmdp and xppimms.
* Added command-line switch to select the source in the ROI for xpmdp.
* --sourceID changed to --srcid throughout.
* xpmdp and xppimms fully refactored.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/314
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/312
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/164


*ixpeobssim (10.1.0) - Tue, 11 Aug 2020 07:39:15 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/108
* Added irgen.mktab.py facility to dump all the ingredients for the MDP
  calculation in tabular format (support for csv and xlsx)


*ixpeobssim (10.0.0) - Mon, 10 Aug 2020 16:36:32 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/106
* New version (v4) of the response function released.
* Generic asymptotic pressure for each of the DU now supported at the IRF
  generation stage---this includes the GPD quantum efficiency, the modulation
  factor and the passive conversion. Version 4 use 687 mbar for all the DUs,
  which is consistent with the sensitivity estimated for the Mission Integration
  Review.
* Combined and non-standard IRFs available in the previous iterations are now
  discontinued.
* xpppims and xpmdp modified to loop aver the three DUs, rather than using
  the combined IRFs.
* New DME density scaling (the same used in ixpesim) and pressure now measured
  in mbar (as opposed to atm) throughout.
* Small refactoring of the code handling the Be window contaminants, and
  certified list Be contaminants from manufacturer is now the default.
* Obsolete files removed and massive cleanup of the IRF documentation.
* pairwise() facility moved into the new module utils.misc.py
* Caching mechanism implemented for loading xcom data.
* A bunch of stuff factored out from irf.modf to evt.mdp.
* Weighted average facility added in ``utils.math_``
* Unit tests added.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/294
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/284
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/295
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/296
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/297
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/298
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/275
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/160


*ixpeobssim (9.0.0) - Mon, 10 Aug 2020 16:27:25 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/103
* Dependence on aplply removed (now relying on astopy.visualizing), docs updated.
* xFITSImage class streamlined and refactored, wit all plotting functionalities
  moved into a base class in core.fitsio.
* Subtle bug fix (pixel offset by one) fixed in binned count maps, as well as
  xFITSImage random sampling---unit tests added.
* Mid-size rafactoring of the xStokesCube class, with arrow-related code
  moved out and rationalized.
* Stokes parameters set to zero outside the physical bounds of the underlying
  interpolator for Stokes sky maps and cubes.
* xpbinview added to the pipeline.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/272
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/281
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/166


*ixpeobssim (8.8.3) - Fri, 07 Aug 2020 08:11:56 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/105
* Fix for skyfield 1.26.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/309


*ixpeobssim (8.8.2) - Fri, 07 Aug 2020 07:40:53 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/104
* Fix for skyfield 1.26.


*ixpeobssim (8.8.1) - Sat, 01 Aug 2020 15:20:02 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/102
* Major cleanup of the hist.py module.
* Unit tests improved.
* pytz dependence removed.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/293


*ixpeobssim (8.8.0) - Sat, 01 Aug 2020 13:53:49 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/101
* Bug fix in xpselect when operating with ds9 region files.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/282
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/300


*ixpeobssim (8.7.0) - Wed, 29 Jul 2020 14:15:29 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/99
* Major restructuring of the support for XSPEC spectral model at simulation time.
* Documentation added.
* Unit tests improved.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/301


*ixpeobssim (8.6.2) - Sat, 25 Jul 2020 16:45:59 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/100
* Fixed a regression triggered by skyfield version 1.23+
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/302


*ixpeobssim (8.6.1) - Wed, 03 Jun 2020 15:28:18 +0200*

* Small facility added for packaging the docs in pdf and zipped html formats.


*ixpeobssim (8.6.0) - Thu, 28 May 2020 21:25:08 +0200*

* Merging in pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/95
* New dependence on skyfield added.
* Basic TLE interface for creating a proxy of the baseline IXPE object.
* Parametrization of the SAA and calculation of the SAA epochs added.
* Calculation of Earth occultation added.
* Calculation of the angles to the Sun and the Moon added.
* A few new facilities in the ``utils.time_`` module added to facilitate time
  conversions.
* Calculation of realistic good time intervals implemented, and all relevant
  pieces of code updated to reflect that (including xpobssim command-line options)
* Initial implementation of a simple visibility tool.
* startmet command-line switch changed to startdate throughout, and default
  observation start changed to 2022-04-21.
* startmet option removed from xpmdp and xppimms.
* Docs added.
* Unit tests added.
* A whola lotta improvements to the docs (branch revamp_docs merged).
* Bonus: fix in the time calculation for periodic sources---the start MET was
  previously ignored.
* Issue(s) closed:
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/232
      * https://bitbucket.org/ixpesw/ixpeobssim/issues/288


*ixpeobssim (8.5.0) - Wed, 20 May 2020 10:52:16 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/96
* Partial major rework of the documentation.


*ixpeobssim (8.4.3) - Sun, 10 May 2020 08:22:24 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/94
* Fix in the documentation.


*ixpeobssim (8.4.2) - Sat, 02 May 2020 07:59:17 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/93
* Confusing option in xpphialign fixed.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/285


*ixpeobssim (8.4.1) - Sat, 02 May 2020 06:56:29 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/89
* Removed any reference to xPolarizationMap class in the tutorial section.
* Fixed a minor typo in source_models: chandra region of interest.
* Added some more details on the xStokesSkyMap class.
* Added in the print statement that the xppims and xpmdp are calculating the
  mdp at the 99 % CL.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/262
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/276


*ixpeobssim (8.4.0) - Thu, 23 Apr 2020 09:29:36 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/91
* xInstrumentalBkg and xPowerLawInstrumentalBkg classes added in srcmode.roi
* New instrumental_bkg config file and example.
* New facilities for projecting GPD positions in the sky.
* Minor refactoring and cleanup.
* Documentation updated.
* Unit tests added.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/264
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/278



*ixpeobssim (8.3.1) - Wed, 22 Apr 2020 18:21:26 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/92
* New mechanism to parse package versions, now handling the PEP 440 specs.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/280


*ixpeobssim (8.3.0) - Wed, 22 Apr 2020 16:13:38 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/90
* Allow changing the default 1--12 keV energy range for xpobssim, via the
  --emin and --emax command-line switches.
* Improved error message when the input spectral model is not positive-definite.
* GRS 1915+105 model completely refactored.
* Generic routine to parse energy-filtered tabular data implemented.
* New facilities to compute the integral flux and integral energy flux in
  srcmodel.spectrum
* Significant improvement in core.spline.scale()
* Documentation modified to reflect the changes.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/267
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/279


*ixpeobssim (8.2.2) - Thu, 27 Feb 2020 11:28:33 +0100*

* Extended the MAX_ENERGY value from 11 to 15 keV in the grs1915_105 config file.

*ixpeobssim (8.2.1) - Mon, 10 Feb 2020 08:56:49 +0100*

* Merging in pull requests
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/88
* XCOM cross sections for H and O added.


*ixpeobssim (8.2.0) - Mon, 10 Feb 2020 08:51:55 +0100*

* Merging in pull requests
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/87
* New infrastructure for the simulation of alpha particles.


*ixpeobssim (8.1.3) - Thu, 30 Jan 2020 13:59:21 +0100*

* Obsolete matplotlib parameters removed from the setup.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/270


*ixpeobssim (8.1.2) - Thu, 30 Jan 2020 12:21:09 +0100*

* Merging in pull requests
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/86
* Fix for CircleCI continuous integration.
* Minor fix to the code for parsing package version strings.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/263


*ixpeobssim (8.1.1) - Thu, 30 Jan 2020 11:23:01 +0100*

* Phony tag because I tagged 8.1.0 on a different branch.


*ixpeobssim (8.1.0) - Thu, 30 Jan 2020 11:18:07 +0100*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/85
* Default value for the BACKFILE and CORRFILE in the header of the SPECTRUM
  extension of binned PHA1* files changed from None to ''
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/254


*ixpeobssim (8.0.0) - Fri, 20 Dec 2019 15:40:37 +0100*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/84
* New round (v3) of response functions.
* New MMA effective area estimated from the post-CDR design
* Effect of all the contaminants in the Be window (as per the Materion data
  sheet) added.
* Photoelectons extracted from the window and the GEM now included in the
  effective area.
* Modulation factor updated based on the calibration of the DU FM2.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/259
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/258
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/257
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/200
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/256
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/255


*ixpeobssim (7.1.0) - Thu, 31 Oct 2019 11:05:58 +0100*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/83
* Added facilities to evt.event to calculate the weighted average of the
  polarization degree for a given photon list.
* Added convenience functions to get the polarization models with the right
  signature (for ixpeobssim) from the Stokes sky-maps and cubes.


*ixpeobssim (7.0.0) - Tue, 29 Oct 2019 09:58:38 +0100*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/81
* Bug fix in the function creating an event mask from a ds9 region file.
* Old ``utils.filter_`` module moved to utils.astro.
* xpselet now accepts arbitrary ds9 region files through the --regfile
  command-line switch (and the temporary regindex switch has been removed.)
* Added facility to retrieve the WCS information from an event file.
* MC_X and MC_Y columns added in MONTE_CARLO extension of the event file.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/181
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/234
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/240


*ixpeobssim (6.5.1) - Mon, 28 Oct 2019 09:27:35 +0100*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/82
* Minor refactoring for the matplotlib color wheels.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/251
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/252


*ixpeobssim (6.5.0) - Fri, 25 Oct 2019 06:39:40 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/80
* Fix in how we calculate the approximate radius of the WCS for Stokes sky maps.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/249


*ixpeobssim (6.4.2) - Thu, 24 Oct 2019 15:36:34 +0200*

* MSH 1552 example updated.


*ixpeobssim (6.4.1) - Thu, 24 Oct 2019 12:13:26 +0200*

* Bug fix in the logger.
* xpphialign added to the reference docs page.
* Small tweak to the top-level Makefile.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/250


*ixpeobssim (6.4.0) - Thu, 24 Oct 2019 11:50:36 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/79
* Layout of the repository Changed
    * ixpeobssim/test -> tests
    * doc -> docs
    * sphinx documentation files refactored and cleaned up
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/195


*ixpeobssim (6.3.0) - Thu, 24 Oct 2019 11:36:15 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/74
* xpphiradalign.py renamed to xpphialign and generalized to handle arbitrary
  polarization models.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/248


*ixpeobssim (6.2.3) - Thu, 24 Oct 2019 05:58:50 +0200*

* Some minor cleanup of old branches.


*ixpeobssim (6.2.2) - Thu, 24 Oct 2019 05:48:06 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/76
* Unit test fixed.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/244


*ixpeobssim (6.2.1) - Thu, 24 Oct 2019 05:10:29 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/74
* Small fix for sphinx 2x
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/246


*ixpeobssim (6.2.0) - Wed, 23 Oct 2019 10:51:25 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/73
* Stokes sky cubes can now load layers from FITS map in the Q/U or polarization
  degree/angle space.


*ixpeobssim (6.1.0) - Tue, 22 Oct 2019 14:01:04 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/71
* Bug fix in the energy interpolation for Stokes sky cubes.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/247


*ixpeobssim (6.0.2) - Tue, 22 Oct 2019 11:34:58 +0200*

* Updated release notes.


*ixpeobssim (6.0.1) - Tue, 22 Oct 2019 11:31:48 +0200*

* This time I tagged the master, but I forgot to update the release notes.


*ixpeobssim (6.0.0) - Tue, 22 Oct 2019 11:27:14 +0200*

(I tagged by mistake a development branch---sorry about that.)

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/71
* Major bug fix in xpphiradalign
* A few bug fixes and improvements in the new xStokesSkyMap and
  xStokesSkyCube classes, which should now be ready for prime time.
* Generic plotting support for polarization models in extended sources
  improved.
* Rotation option added when loading Stokes sky maps from polarization
  degree and angle.
* Cleanup of the config files for extended sources, and files renamed.
* Polarization map arrows now working on plain matplotlib figures.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/241
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/188
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/194
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/245
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/235
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/243


*ixpeobssim (5.4.0) - Fri, 18 Oct 2019 16:44:37 +0200*

(Mind this was tagged for logistical reasons with known issues and you are
strongly advised against using this version. Please update straight to 6.0.0)

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/70
* New Tycho config file and example using the Stokes sky-cube functionality.
* Improvements to xpphiradalign: new command-line options, and fix in the
  pipeline wrapper and alignment algorithm refactored.
* Error message prompted upon sky coordinates outside the Stokes sky-map grids
* Logger tweaked to print messages in color.


*ixpeobssim (5.3.0) - Wed, 16 Oct 2019 09:35:10 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/69
* New small tool to align the photoelectron directions radially
* Tool made available into the pipeline.


*ixpeobssim (5.2.1) - Wed, 16 Oct 2019 09:18:45 +0200*

* Custom argument formatter added to the argument parser to allow for more
  flexibility in the parser output.
* __description__ for all the applications revamped.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/239


*ixpeobssim (5.2.0) - Tue, 15 Oct 2019 11:02:39 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/67
* Added an option --innerrad to select annuli in xpselect
* Added an option --invert to invert a selection in xpselect
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/237
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/238


*ixpeobssim (5.1.0) - Mon, 14 Oct 2019 12:31:37 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/67
* Full refactor of all the classes dealing with polarization maps into a single
  xStokesSkyMap class, that is able to read in FITS files with Q/U, x/y, or
  polarization degree/angle, but only uses the Stokes parameters internally.
* Configuration files and examples updated.
* New xStokesSkyCube added, for 3-d interpolation in x, y and energy.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/119
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/181


*ixpeobssim (5.0.0) - Wed, 09 Oct 2019 12:26:27 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/66
* GEM charging introduced in the simulation.
* GEM gain added to the output file, in the MONTE_CARLO extension.
* xpchrgcorr.py added to process event lists and correct for the charging.
* New EFLUX binning algorithm implemented.
* Pulse-height analysis code refactored (issue #224)
* Simple decorator to time functions added.
* Horrible inverse spline for the energy to channel conversion replaced with
  a binary search.
* Energy dispersion classes modified to keep track of the measured energy
  before digitization
* Various minor refactorings and fixes.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/231
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/221
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/224
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/71


*ixpeobssim (4.2.1) - Tue, 01 Oct 2019 10:52:55 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/65
* Minor update to the data format (some D changed to E)
* Definition of the event-by-event Stokes parameters reverted back to the
  original definition (with the extra factor of 2).
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/218
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/220


*ixpeobssim (4.2.0) - Wed, 18 Sep 2019 09:45:51 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/64
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/219


*ixpeobssim (4.1.0) - Mon, 19 Aug 2019 17:39:40 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/63
* Initial implementation of a spiral dithering pattern (unit-tested but not
  yet used for either analysis or simulations).


*ixpeobssim (4.0.0) - Thu, 04 Jul 2019 16:01:47 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/62
* Dithering pattern implemented in the simulation and fully exposed through
  ixpeobssim command-line switches.
* DETPHI output column is now correctly rotated with respect to PHI on a
  DU by DU basis.
* GPD fiducial area is now properly applied, and the DU rotation is implemented
  in conformance with the relevant IXPE mechanical interface control document.
* Unit tests added.
* Typo in a function name in utils.units fixed.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/193


*ixpeobssim (3.6.0) - Fri, 28 Jun 2019 13:32:51 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/60
* This is mainly to allow disengaging the application of the vignetting for
  situations where the source is generated into the detector and the optics
  are irrelevant (e.g., the internal backgroud). But really, this is a
  significant refactoring of the entire srcmodel.roi module.
* Added a base class for model components representing calibration flat
  fields.
* Added check to make sure that the event times are sorted at the event-list
  filling time.
* xEventList.__len__() deprecated.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/185
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/214
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/215


*ixpeobssim (3.5.1) - Thu, 27 Jun 2019 11:25:15 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/61
* Patch to avoid tracebacks from xFITSImage
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/167
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/213


*ixpeobssim (3.5.0) - Wed, 26 Jun 2019 12:52:22 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/58
* Added a check that the pdf is positive for all the random number generators.
* Added diagnostics for unphysical polarization values.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/39
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/174


*ixpeobssim (3.4.0) - Wed, 26 Jun 2019 12:20:24 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/57
* This is essentially an update of the data format, to align it more with the
  corresponding document.
* Unrelated minor fix in utils.system.cmd(), see issue #211
* Issue closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/175
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/206
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/207
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/211


*ixpeobssim (3.3.0) - Wed, 26 Jun 2019 12:17:55 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/59
* Random seed set to a random value by default.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/189


*ixpeobssim (3.2.0) - Sun, 23 Jun 2019 08:54:24 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/55
* Refactoring of the spectrum model, where the time-dependent spectral
  parameters are now handled in a sensible and consistent fashion.
* Issue closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/202


*ixpeobssim (3.1.0) - Wed, 19 Jun 2019 14:24:20 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/45
* Added a tutorial section to the documentation.


*ixpeobssim (3.0.0) - Wed, 19 Jun 2019 14:17:16 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/52
* Chandra2ixpe application removed to avoid potential confusion between xpobssim
  and the converter (the two applications were almost identical). As a result of
  this change, a standard simulation or a conversion from Chandra will take
  place based on the type of ROI implemented in the source configuration file.
* Issue closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/192


*ixpeobssim (2.10.0) - Wed, 19 Jun 2019 12:37:11 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/56
* xpeInterstellarAbsorptionModel changed to xInterstellarAbsorptionModel
* README file updated, and all contents from the wiki removed.
* Programmatic import cleaned up to fix a warning with Python 3.
* Copyright notice in the startup message updated.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/203
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/70
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/165
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/198


*ixpeobssim (2.9.0) - Tue, 18 Jun 2019 16:25:04 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/54
* This is a massive refactoring of all the spline classes, including those
  dealing with the generation of random numbers, the two most profound changes
  being a fix of the way we handle meshgrids in bivariate splines, and a change
  in the signature of random number generators with auxiliary variables
  (rv and aux swapped). Both changes involve a subtle variation of the semantics
  but the net result is now clean and self-consistent, and the operation was
  absolutely necessary before we started implementing new complex features on
  top of our spline/rng ecosystem. (This should be all largely transparent to
  end user but, needless to say, we might have introduced side effects, so
  watch out.)
* xSourceSpectrum class added.
* As a consequence of the spline refactoring our interbale representations of
  the energy dispersion and the azimuthal response generators, as well as the
  source and count spectra, are now rotated by 90 degrees.
* Interpolation in log space now working.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/196
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/163
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/88
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/54
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/58
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/113
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/126
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/115
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/197
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/199
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/201


*ixpeobssim (2.8.2) - Mon, 03 Jun 2019 14:45:08 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/50
* Two minor changes to the data format to conform to the FITS standards
  ("DEC---TAN" changed to "DEC--TAN" and comment for the "DATE" field in the
  event lists shortened).


*ixpeobssim (2.8.1) - Mon, 03 Jun 2019 14:39:44 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/49
* One unit test on the unpolarized response added (this was triggered by
  issue #182.)


*ixpeobssim (2.8.0) - Mon, 03 Jun 2019 14:21:44 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/48
* "Bug" in the time assignment when converting Chandra observations (see
  issue #191---up to now we were using the Chandra event times and this, in
  turn,was causing issues with the CCD data). The Chandra event times are now
  ignored, and assigned randomly between the start and the end of the
  observation, instead.
* Unit tests added.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/191


*ixpeobssim (2.7.0) - Mon, 03 Jun 2019 10:56:08 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/48
* Top-level Makefile changed, now it works on different devices
* New Furier_harmonics method to build pulse profiles and/or
  phase-dependent polarization degree function in the configuration
* New cutoff_power_low in spectrum.py
* Implemented xXspecModel and xXspecModelBuilder classes for a generic
  source configuration, including PYXSPEC support


*ixpeobssim (2.6.1) - Tue, 28 May 2019 18:20:32 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/46
* This fixes an issue with the trigger ID column for periodic sources.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/142


*ixpeobssim (2.6.0) - Tue, 02 Apr 2019 14:57:50 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/41
* Major revision of the distutil-related files---now we should be up and
  and running into user-mode.
* Issue(s) closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/176


*ixpeobssim (2.5.0) - Wed, 07 Nov 2018 07:03:51 +0100*

* Merged in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/43
* Merged the polarization class with the combined polarization class.
* Removed the xpolarizationarrows class.
* Added a draw arrows method to the xBinnedStokesCube class.
* Added the vela nebula pulsar configuration and example script
* Fixed the msh1552 configuration and example script.


*ixpeobssim (2.4.0) - Wed, 10 Oct 2018 09:09:23 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/40
* XFLT0001 keyword added to the SPECTRUM extension of the binned Stokes spectra.
* Added support for native spectro-polarimetric fit in XSPEC.


*ixpeobssim (2.3.0) - Sun, 16 Sep 2018 22:31:51 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/37
* Stokes-spectra output now implemented in xpbin.
* New .mrf response functions in the caldb.
* evt.xspec module update to support spectro-polarimetric fitting in XSPEC.
* Option mc=True for xpbin.py in PHA1* mode correctly implemented.
* Pipeline now saving and plotting figures automatically.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/155,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/172


*ixpeobssim (2.2.1) - Thu, 13 Sep 2018 16:50:39 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/39
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/173


*ixpeobssim (2.2.0) - Wed, 12 Sep 2018 21:55:45 +0200*

* Merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/38
* Added facilities to profile the memory usage.
* Avoid copying HDU lists when selecting rows within an event file, which in
  turn fixes a serious memory leak in xpselect, see https://bitbucket.org/ixpesw/ixpeobssim/commits/6a7b62cc55856892f751c55b7361b9531b9027a6
* Pipeline now saving and plotting figures automatically.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/149,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/172


*ixpeobssim (2.1.1) - Mon, 03 Sep 2018 16:53:22 +0200*

* Minor changes, tagging the correct branch (master).


*ixpeobssim (2.1.0) - Mon, 03 Sep 2018 16:42:47 +0200*

* Full rewrite of the pipeline module, which is now a plain set of function,
  as oposed to a class.
* Documentation of the new pipeline module added.
* All examples reviewed in light of the new pipeline framework.
* Application signatures in ixpeobssim/bin uniformed to make it easier to wrap
  them.
* pairwise() and pairwise_enum() methods added to binning.
* POL_ANGLE and POL_ANGLE_ERR fields for binned modulation cubes changed to
  POL_ANG and POL_ANG_ERR for uniformity with the rest of the framework.
* Default DU suffix for the output file name changed from 1, 2, 3, to
  du1, du2, du3.
* Full rewrite of the old evt.fitting module, which is now evt.xspec; like
  the new pipeline, the module is now a collection of functions, as opposed to
  a class.
* Major revamp of the crab_pulsar configuration file and associated example
  pipeline, with docs added in the proper section.
* Major revamp of the J1708 configuration file and associated example
  pipeline, with docs added in the proper section.


*ixpeobssim (2.0.0) - Sun, 26 Aug 2018 15:25:59 +0200*

* Version v2 of the IXPE response functions released, with many changes, see https://bitbucket.org/ixpesw/ixpeobssim/issues/161
      * IRFs are now defined from 1 to 12 keV (was 15 in the previous
        iteration).
      * new MMA effective area parametrization, with thicker thermal shields,
        see https://bitbucket.org/ixpesw/ixpeobssim/issues/152
      * effect of the aluminination of the Be window included, see
        https://bitbucket.org/ixpesw/ixpeobssim/issues/151
      * transparency of the UV filter added, see
        https://bitbucket.org/ixpesw/ixpeobssim/issues/159
      * GPD quantum efficiency updated with the tabulated DME density at 20
        degrees C.
* Added new modulation response function in the CALDB, and associated tools to
  load and visualize it.
* Added toy response functions for debugging purposes.
* Major refactoring of the code generating the response functions:
     * xcom interface added, and all the relevant GPD characteristics are
       now calculated on the fly, rather than written in text files.
     * calibgen folser moved to irfgen
     * all irfgen data files reorganized into folders matching the
       new irfgen Python modules.
* Major refactoring of the IRF classes:
     * xIrfBase class created, and now all the response-file classes are
       inheriting from it.
* Major refactoring of the command-line parsing facilities.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/161,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/152,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/153,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/151,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/159


*ixpeobssim (1.8.1) - Fri, 24 Aug 2018 15:25:48 +0200*

* ngc1068 configuration file (and pipeline) adapted to the new source
  etiquette.
* Source documentation tweaked.


*ixpeobssim (1.8.0) - Wed, 22 Aug 2018 16:55:03 +0200*

* Merged in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/
* All toy models renamed and revamped following a minimal set of guidelines.
* Pipelines added for all toy models, with the corresponding documentation
  linked from the main docs page.
* Section about source etiquette added in the documentation.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/141


*ixpeobssim (1.7.0) - Fri, 17 Aug 2018 14:07:57 +0200*

* Merged in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/32
* Active area of the GPD readout on the focal plane implemented (events
  outside the detector are trimmed out).
* DU rotation implemented (each of the three DU now gets its own orientation).
* Added option to set the telescope roll angle in xpobssim and chandra2ixpe.
* binsz (now pixsize) functionality restored in xpbin.
* Default image size for xpbin in CMAP mode adjusted.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/48


*ixpeobssim (1.6.0) - Wed, 08 Aug 2018 15:10:03 +0200*

* Major refactoring of the binning classes, merging in pull reques
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/31
* Sum of the relevant binning data structures now implemented for all types
  via the __iadd__ (as opposed to __add__) overload.
* PHASG binning algorithm renamed as PP.
* nxpix and nypix command-line switches for binned map now unified in npix
  (i.e., we shall be only producing square maps).
* x-axis errors supported in xScatterPlot.
* Modulation cube objects now using the internal information for
  visualization purposes (as opposed to fitting).
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/93,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/120,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/130,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/134,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/145


*ixpeobssim (1.5.1) - Thu, 02 Aug 2018 15:48:13 +0200*

* Cleanup imports, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/29
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/140


*ixpeobssim (1.5.0) - Wed, 01 Aug 2018 21:24:22 +0200*

* Windows support added, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/29


*ixpeobssim (1.4.1) - Wed, 01 Aug 2018 21:05:44 +0200*

* merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/22
  (this only involves files in the sandbox)


*ixpeobssim (1.4.0) - Wed, 01 Aug 2018 21:03:13 +0200*

* Improved azimuthal response generator, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/25
* This is closing issues https://bitbucket.org/ixpesw/ixpeobssim/issues/49 and
  https://bitbucket.org/ixpesw/ixpeobssim/issues/135


*ixpeobssim (1.3.0) - Wed, 01 Aug 2018 20:46:17 +0200*

* NGC 1068 model added, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/28
* Full implementation of the harmonic addition theorem added in
  srcmodel/polarization.py


*ixpeobssim (1.2.0) - Wed, 01 Aug 2018 20:38:54 +0200*

* Added entry points for the user installation, merged in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/24
* Default output folder changed to ~/ixpeobssimdata
* xpbinviewer.py and xpirfviewer.py renamed to xpbinview.py and xpirfview.py
* Some obsolete scripts in ixpeobssim/bin moved to the sandbox
* Documentation updated.


*ixpeobssim (1.1.0) - Thu, 26 Jul 2018 12:55:32 +0200*

* Tycho example added, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/26
* Emergency fix for issue https://bitbucket.org/ixpesw/ixpeobssim/issues/137


*ixpeobssim (1.0.0) - Tue, 17 Jul 2018 15:26:36 +0200*

* New set of irfs available for each detector unit with 350 energy channels.
* Added support to the new caldb structure.
* All applications, test scripts and examples updated to support the new
  du-based irfs.
* Maps of polarization degree and angle can now be used to build the input
  polarization model.
* Added methods to plot and overlay arrows to the input and output counts maps.
* New ARMAP binning algorithm implemented.
* Major refactoring of Stokes cubes data format.
* add() and from_list() methods implemented in (almost) all binning classes.
* Added a flag to exclude regions in xChandraObservation class.
* Added a simulation of msh1552 using models provided by N. Bucciantini.
* Added a simulation of AXP 1RXS_J1708 using models provided by R. Turolla and
  R. Taverna.
* Merged in pull request:
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/23
* Added a new application to calculate the phase and append the column to the
  event list. Merging in pull request:
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/19
* Phase column removed from output file of xpobssim and chandra2ixpe (issue #9).
* Fixed crab_pulsar, casa, cena examples. Merging pull requests:
    * https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/15,
    * https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/18,
    * https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/21.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/66,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/67,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/68,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/80,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/121,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/125,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/116.
* Added sandbox folder.


*ixpeobssim (0.64.0) - Mon, 21 May 2018 12:27:08 +0200*

* Bug fix in the simulation of periodic sources, marging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/16
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/118


*ixpeobssim (0.63.1) - Fri, 27 Apr 2018 15:03:38 +0200*

* __future__ imports added everywhere to support Python 2, merging in pull
  request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/14/


*ixpeobssim (0.63.0) - Sat, 14 Apr 2018 09:28:20 +0200*

* Creator ID added to all the output fits files.
* chandra2ixpe unit test and pipeline interface fixed.
* Merged in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/13
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/104


*ixpeobssim (0.62.1) - Sat, 14 Apr 2018 01:47:17 +0200*

* Visibility changed to modulation throughout, and some docs added, merged in
  pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/12/fix-visibility/diff
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/109


*ixpeobssim (0.62.0) - Fri, 13 Apr 2018 20:03:01 +0200*

* Data format updated to the new, tentative definition of the Level 2 files,
  following the I2T face to face meeting on April 13, 2018, mergin in pull
  request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/11
* Highlights include the fact that output photon lists can now be
  correctly parsed by the standard tools of the community, including
  ds9 and ximage.
* Deadtime-related keywords included in all the relevant headers, and
  propagated to the PHA1 binned files for spectral fitting with xspec.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/92,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/8


*ixpeobssim (0.61.0) - Wed, 11 Apr 2018 11:49:33 +0200*

* Major refactoring of the modeling, fitting and plotting classes, merging in
  pull request https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/10
* Major restructuring of the modulation cube class.
* Refactoring of the binning module.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/34,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/56,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/105,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/106


*ixpeobssim (0.60.0) - Sat, 07 Apr 2018 10:53:02 +0200*

* xpselect refactoring for the new data format, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/9
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/102
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/103


*ixpeobssim (0.59.0) - Fri, 06 Apr 2018 14:01:24 +0200*

* Complete refactoring of the chandra2ixpe converter, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/8
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/20,
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/81
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/85
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/89
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/91


*ixpeobssim (0.58.0) - Tue, 27 Mar 2018 23:49:02 +0200*

* Tentative implementation of the new data format and major refactoring
  of the xpobssim internals, merging pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/7
* More unit tests for the simulation and analysis pipelines.
* Issues closed:
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/90
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/94
    * https://bitbucket.org/ixpesw/ixpeobssim/issues/38


*ixpeobssim (0.57.0) - Sat, 24 Mar 2018 10:24:07 +0100*

* Documentation revamped and improved, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/5
* A typo in all the header files fixe, see issue #5, merging in pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/6


*ixpeobssim (0.56.0) - Sun, 11 Mar 2018 18:12:06 +0100*

* Now officially supporting Python 3x (and discontinuing Python 2x). Merged in
  pull request
  https://bitbucket.org/ixpesw/ixpeobssim/pull-requests/4


*ixpeobssim (0.55.1) - Mon, 02 Oct 2017 16:53:25 +0200*

* Minor.


*ixpeobssim (0.55.0) - Mon, 02 Oct 2017 16:15:34 +0200*

* A few new models added or updated (Crab, Cas A, MSH 15-52).
* Added code to generate IRFs in the new CALB format (not yet used).


*ixpeobssim (0.54.1) - Tue, 01 Aug 2017 08:53:52 +0200*

* More work on the documentation.


*ixpeobssim (0.54.0) - Mon, 31 Jul 2017 10:11:13 +0200*

* Extensive restructuring of the documentation.


*ixpeobssim (0.53.0) - Sat, 29 Jul 2017 09:25:35 +0200*

* Added an example script (examples/uniform_disk_stokes.py) to plot the
  polarization degree and angle using the stokes cube.
* Added a script (ximpol/detector/gpd.py) and associated files to cross check
  the calculation of the window transmission and the GPD efficiency.
* Added SCUBE algorithm to xpbin.
* Added new events binning classes to make and read the stokes cubes (SCUBE).
* Added method to use stokes parameters to build the polarization maps in
  example script casa_pol_map.py.
* Added xppimms to the pipeline.
* A couple of modifications introduced by mistake reverted.
* xppimms and chandra2ximpol scripts made executable.


*ximpol (0.52.0) - Tue, 24 Jan 2017 12:14:17 +0100*

* First implementation of xppimms for MDP calculation with source paramteres
  provided through command-line switches (issue #153)
* Added example script to compare polarization values found via mcube vs stokes,
  using a single point source.
* Added polarization fraction method to class xStokesAccumulator
* Refactoring of evt.fitting.xSpectralFitter class.
* Implemented the vignetting in the chandra2ximpol converter.
* Added a function to draw the psf circle in irf.psf.xPointSpreadFunction
  (closing issue #151).
* Scripts to simulate GRS1519+105 and the txt files provided by Banafsheh.
* Minor refactoring of the IRF plotting code (closing issue #150).
* Some tweaks to the xpmdp.py script (closing issue #107).
* Bug fix for issues #25 and #22 (closing them).
* Observation plan plot added in the examples folder.
* Scripts to simulate GRS1519+105 using the txt files provided by Michal.
* Added analysis pipeline to simulate Cen A using chandra-to-ximpol converter.
* Added a first implementation of a "Stokes accumulator" class (see issue #148).


*ximpol (0.51.1) - Fri, 01 Jul 2016 14:04:20 +0200*

* Bug fix in the check for the polarization degree (see issue #73).


*ximpol (0.51.0) - Fri, 01 Jul 2016 12:00:53 +0200*

* First stub at implementing splines in log space (issue #20).
* Complete PSF refactoring.
* New sets of IRFs created, and default now pointing to
  'xipe_mirror-30s-f4_psf-jetx-rescaled-hew30_gpd-baseline'


*ximpol (0.50.0) - Thu, 30 Jun 2016 14:56:31 +0200*

* Blazar sensitivity plot added in the examples folder.
* New config file and associated pipeline added for J1708.
* Avoid reading the modulation factor in xBinnedModulationCube
  (using the effective factor written in the fits file instead, when it
  comes to converting a visibility into a polarization fraction.)
* One unit test added.
* New bivariate spline class added supporting orders greater than 1 on both
  axes.
* Bug fix for issue #143 (closing it).
* xUnivariateAuxGenerator now supporting spline orders greater than 1 on
  both axes (taking advantage of the new bivariate spline class).
* Fix for issue #73 (clsing it).


*ximpol (0.49.0) - Fri, 24 Jun 2016 16:31:55 +0200*

* "xipe_goal" IRFs now used by default by all the xp tools.


*ximpol (0.48.0) - Wed, 22 Jun 2016 23:12:31 +0200*

* Added script to make the polarization map of casa in examples
* Added the option to draw the psf circle in the count map of casa in the main
  casa.py example.
* Added the Cyg-X1 ascii files for the model at 40 degree inclination.
* Added checks in the univariate spline constructors to make sure that the
  input x-values are sorted and unique (closing issue #84.)
* binning module modified in such a way that all the relevant quantities
  are also calculated over the entire energy range when the energy
  binning is longer than 1 (i.e., if you do 2--4 keV and 4--8 keV, you
  get 2--8 keV for free).


*ximpol (0.47.3) - Wed, 22 Jun 2016 08:56:35 +0200*

* Added M87 configuration files and analysis pipeline for chandra-to-ximpol
  converter.
* Added configuration files and analysis pipeline example for MCG-6-30-15,
  ARK 120 and NGC 1365.
* Configuration file for GK Per slightly tweaked.


*ximpol (0.47.2) - Sat, 18 Jun 2016 06:56:49 +0200*

* Minor doc update.
* Added GK Per configuration file and analysis pipeline example.


*ximpol (0.47.1) - Thu, 16 Jun 2016 16:39:40 +0200*

* Slight tweak to the command-line switches for xpbin and xpselect, in order
  to be able to set the mc keyword argument to True from the pipeline.
* Fix for issue #105 (closing it).


*ximpol (0.47.0) - Thu, 16 Jun 2016 15:34:41 +0200*

* Redshift added to the simulation (issue #121).


*ximpol (0.46.0) - Thu, 16 Jun 2016 15:26:42 +0200*

* Interstellar absorption added to the simulation (closing issue #120).


*ximpol (0.45.1) - Wed, 15 Jun 2016 15:00:50 +0200*

* Minor bug fix.
* ximpol.srcmodel.polarization module added to the documentation.


*ximpol (0.45.0) - Wed, 15 Jun 2016 14:46:01 +0200*

* Galactic absorption column density and redhsift added to all the model
  components (working our way through issue #120).
* Some cleanup of the Galactic absorption code (issue #120).
* Unit test for the intertellar absorption added.
* Fix for issue #103.
* Some significant refactoring of the xpmdp code.
* xpmdp added to the pipeline (closing issue #137).
* Added a unit test comparing the xpmdp.py output with the figures from the
  online sensitivity calculator.


*ximpol (0.44.2) - Fri, 10 Jun 2016 17:58:27 +0200*

* Trivial fix.


*ximpol (0.44.1) - Fri, 10 Jun 2016 17:54:30 +0200*

* Installation notes updated (closing issue #131).


*ximpol (0.44.0) - Fri, 10 Jun 2016 15:35:16 +0200*

* Added method to the xPolarizationMap object to plot the polarization
  degree and angle (in addition to the x and y components).
* Added new corona models for Cyg-X1 in ascii folder
* Added new Cyg-X1 config files (one per model) in the config folder
* Added new script tp make the plot comparing the polarization fraction for
  Cyg-X1 for the two coronal models in examples
* Added method to run several simulations (with different seeds) and merge the
  output files to one single. This is in the run method for the Cyg-X1 example.
* Added script to make the map of the polarization degree and map of the sigma
  for the polarization degree.
* Initial import of the module and files for parametrizing the Galactic
  absorption (issue #120).
* Quick fix to issue #129.


*ximpol (0.43.0) - Sat, 04 Jun 2016 22:13:08 +0200*

* Added a script to plot MDP for the crab pulsar with the nebula background in
  examples.
* Added effective mu, source counts and mdp to the xBinnedModulationCube class
* Showcase updated with the fix to the PSF.
* First implementation of the chandra2ximpol converter.


*ximpol (0.42.1) - Wed, 20 Apr 2016 21:11:25 +0200*

* Team updated.


*ximpol (0.42.0) - Wed, 20 Apr 2016 11:30:41 +0200*

* xpmdp.py adapted (via brute-force) to periodic sources.
* MDP information added to the xpbin output in MCUBE mode.


*ximpol (0.41.1) - Tue, 12 Apr 2016 16:54:48 +0200*

* Minor fix in the MDP calculator output.


*ximpol (0.41.0) - Tue, 12 Apr 2016 16:28:11 +0200*

* New utils/units module added.
* Configuration file for Abell 85 (along with its configuration file) added.
* Source model string formatting improved (issue #101).


*ximpol (0.40.0) - Tue, 12 Apr 2016 13:38:08 +0200*

* Initial implementation of the Cyg-X1 config/example.
* First implementation of a script for the calculation of the MDP.
* Significant refactoring of the srcmodel/spectrum.py module.


*ximpol (0.39.4) - Thu, 07 Apr 2016 13:59:05 +0200*

* Cas A example tweaked.


*ximpol (0.39.3) - Thu, 07 Apr 2016 07:00:47 +0200*

* Added minimal support for log scale on the the z axis when plotting
  bivariate splines.
* Added doc/scripts folder (work in progress).


*ximpol (0.39.2) - Tue, 22 Mar 2016 14:51:31 -0700*

* Cas A movie updated.


*ximpol (0.39.1) - Tue, 22 Mar 2016 13:59:13 -0700*

* Added a GRB example to the gallery.
* Added the Cas A movie.


*ximpol (0.39.0) - Fri, 18 Mar 2016 11:54:43 -0700*

* GRB 130427 configuration file revamped.
* Bug fix in the binning module for the LOG time binning.
* Some more infrastructure in place for arbitrary source-based sampling times
  (issue #44).
* Added a new example for GRB 130427.


*ximpol (0.38.1) - Thu, 17 Mar 2016 09:51:51 -0700*

* References for the Crab pulsar example added.


*ximpol (0.38.0) - Wed, 16 Mar 2016 15:48:56 -0700*

* One more unit test added.
* A few tweaks and some cleanup.
* Optional scale and offset parameters added to the plot() method for the
  univariates splines.
* Bug fix for issue #97.
* Crab pulsar example revamped.


*ximpol (0.37.1) - Tue, 15 Mar 2016 17:10:54 -0700*

* Crab pulsar added to the showcase.


*ximpol (0.37.0) - Tue, 15 Mar 2016 15:15:36 -0700*

* Added a pipeline example for the Crab pulsar.
* Equipopulated-binning code refactored (issue #93).
* evt.select.py renamed as evt.subselect.py (issue #96).
* xpxspec refactored, with most of the code being moved to evt.fitting.py
  (issue #92).
* Some specific refactoring.
* Equipopulated binning refactored (issue #93).


*ximpol (0.36.1) - Sat, 12 Mar 2016 07:40:17 -0800*

* First complete Cas A section in the gallery (issue #80).


*ximpol (0.36.0) - Sat, 12 Mar 2016 06:03:59 -0800*

* Initial stub at the ximpol gallery (issue #80).
* Short version of the command-line switches removed from xpxspec, and all of
  them passed as keyword arguments (issue #71).
* xpxspec added to the pipeline, and a new example added.
* More tweaks to the Cas A analysis pipeline example.
* ebinning LIST mode added to xpbin.
* Significant refactoring of the xBinnedModulationCube
  class to allow to reuse single analysis/plotting tasks externally.
* Pretty much done with the lamp_post pipeline example.
* A few interface tweaks.
* Fix for issue #77.
* Getting started on documenting the architecture of the package.
* Model for Tycho added.
* Bug in the PSF fixed (issue #82).
* A few files renamed (removed the leading test) to prevent issues with the
  unit testing.


*ximpol (0.35.4) - Wed, 09 Mar 2016 16:54:43 -0800*

* Enforce data-type consistency in the output event files (issue #66).


*ximpol (0.35.3) - Wed, 09 Mar 2016 14:01:52 -0800*

* Clobber mechanism implemented at the single tool level and propagated to
  the pipeline.


*ximpol (0.35.2) - Wed, 09 Mar 2016 10:06:32 -0800*

* Smoother version of the Cas A spectral models.


*ximpol (0.35.1) - Tue, 08 Mar 2016 22:01:46 -0800*

* New Cas A spectral models.


*ximpol (0.35.0) - Tue, 08 Mar 2016 15:31:58 -0800*

* Modified configuration file for Cas A, now with separate extended
  components for the thermal and the non-thermal emission.
* Subtraction implemented for unidimentional splines.
* Classes for the fit of the azimuthal distributions tweaked.
* One full analysis pipeline for Cas A implemented in ximpol/examples/casa.py.
* A few obsolete files removed.


*ximpol (0.34.1) - Mon, 07 Mar 2016 16:54:42 -0800*

* Help formatter for xpobssim, xpselect and xpbin changed.


*ximpol (0.34.0) - Mon, 07 Mar 2016 16:45:19 -0800*

* Internals of xpobssim, xpselect and xpbin tweaked to be fully configurable
  via keyword arguments, to facilitate pipelining analyses.
* ximpol/examples folder added.
* First stup at an analysis pipeline to facilitate complex simulation/analysis
  chains.


*ximpol (0.33.1) - Mon, 07 Mar 2016 12:11:25 -0800*

* Minor changes to the lamp_post_accreting_bh source model.


*ximpol (0.33.0) - Mon, 07 Mar 2016 10:59:30 -0800*

* Some changes for the creation of the Cas A movie.
* Minor tweaks.
* Layout of the configuration files reorganized, with ximpol/srcmodel/config
  moved to ximpol and the ascii and fits folder moved as subfolders therein.
* ximpol.__init__.py and cleanup script modified accordingly.
* Lamp-Post accreting BH model (by Alberto) added.
* All soure models adapted to the new layout.
* And all the unit tests run for cross-check.


*ximpol (0.32.1) - Wed, 02 Mar 2016 10:44:26 -0800*

* Changed name of xpbinview.py to xpviewbin.py.
* Short options removed from xpobssim (issue #71).


*ximpol (0.32.0) - Tue, 01 Mar 2016 16:45:31 -0800*

* Implemented position-dependent polarization patterns based on FITS maps.
* All configuration files updated to the new interfaces.


*ximpol (0.31.0) - Tue, 01 Mar 2016 15:17:37 -0800*

* Replace numpy.fill() with numpy.full() when appropriate (issue #66).
* New display interface to binned files ximpol/bin/xpview.py (issue #55).
* Obsolete script ximpol/bin/xpimgview.py removed.
* Obsolete script ximpol/bin/xpevtview.py removed.
* A couple of bug fixes in the source models.


*ximpol (0.30.1) - Sat, 27 Feb 2016 08:20:34 -0800*

* Closed issue #63.


*ximpol (0.30.0) - Sat, 27 Feb 2016 06:55:48 -0800*

* A couple of command-line switches added to xpselect (issue #51).
* xpbin options propagated to the output files (issue #60).


*ximpol (0.29.0) - Fri, 26 Feb 2016 18:41:42 -0800*

* Source model for Cas A added.
* First xpselect implementation (issue #51).
* Subtle bug fix in the CMAP binning (issue #70).


*ximpol (0.28.1) - Thu, 25 Feb 2016 16:48:43 -0800*

* Updated installation instructions (issue #64).


*ximpol (0.28.0) - Thu, 25 Feb 2016 15:51:26 -0800*

* Phaseograms implemented in xpbin.py (issue #67).


*ximpol (0.27.0) - Thu, 25 Feb 2016 15:31:53 -0800*

* Work started toward the implementation of periodic sources (issue #43).
* New xEphemeris class in ximpol.srcmodel.roi.py
* New xPeriodicPointSource class in ximpol.srcmodel.roi.py
* Some significant refactoring of the spline and rand classes to allow for
  more flexibility.
* Major change to the source model interface---the energy spectrum and
  polarization degree and angle are now passed to the constructor.
* A whole bunch of obsolete stuff removed from ximpol.srcmodel.spectrum
  (issue #64).
* All configuration files reworked according to the new interfaces.


*ximpol (0.26.0) - Tue, 23 Feb 2016 16:42:27 -0800*

* FILE mode implemented for tbinalg (issue #53).


*ximpol (0.25.0) - Tue, 23 Feb 2016 16:33:27 -0800*

* ebinalgs FILE and EQP implemented (issue #56).


*ximpol (0.24.1) - Tue, 23 Feb 2016 15:55:06 -0800*

* Fixed unit tests.


*ximpol (0.24.0) - Fri, 19 Feb 2016 16:14:36 -0800*

* Vignetting now into the effective area tables (but not used in the
  simulation, yet).


*ximpol (0.23.1) - Thu, 18 Feb 2016 15:03:59 -0800*

* More information added to the IRF primary headers (issue #49).


*ximpol (0.23.0) - Thu, 18 Feb 2016 14:56:15 -0800*

* Major refactoring of ximpol/detector/xipe.py to use the new classes
  (issue #49).
* New optics aeff files provided by Fabio committed (but only the on-axis
  values used for the time being).
* XIPE baseline and goal response functions created (only the effective areas
  differ for the time being).


*ximpol (0.22.4) - Mon, 08 Feb 2016 16:34:11 -0800*

* Fix for issue #59.


*ximpol (0.22.3) - Mon, 08 Feb 2016 16:25:59 -0800*

* Fix for issue #58.


*ximpol (0.22.2) - Mon, 08 Feb 2016 15:51:53 -0800*

* Quick polarization analysis routine in place.
* Bug fix in the new code reading the IRFs.


*ximpol (0.22.1) - Mon, 08 Feb 2016 15:11:38 -0800*

* More refactoring of the binning classes.
* Detector, ROI and IR information propagated from the event to the binned
  files (issue #57).


*ximpol (0.22.0) - Fri, 05 Feb 2016 13:56:10 -0800*

* MCUBE mode implemented in xpbin.py


*ximpol (0.21.2) - Thu, 04 Feb 2016 15:41:41 -0800*

* Source model string formatting improved.
* A few minor changes.


*ximpol (0.21.1) - Thu, 04 Feb 2016 14:28:43 -0800*

* Committed a whole bunch of files left out by mistake.


*ximpol (0.21.0) - Thu, 04 Feb 2016 14:27:20 -0800*

* Major refactoring and revamp of xpevtview.py
* New class for tabulated stationary spectra.
* New configuration file for the SgrB complex.
* Spectral data for the SgrA and SgrB complexes.
* New small utility (xpsrccoords.py) to search for source coordinates.


*ximpol (0.20.0) - Thu, 04 Feb 2016 10:43:26 -0800*

* Gaussian disk spatial template implemented.
* A few srcmodel config files renamed.


*ximpol (0.19.1) - Wed, 03 Feb 2016 16:17:09 -0800*

* Updated documentation.


*ximpol (0.19.0) - Wed, 03 Feb 2016 16:12:42 -0800*

* Uniform disk implemented (issue #54).
* Added command-line option to use the MC Ra/Dec for xpbin.


*ximpol (0.18.0) - Wed, 03 Feb 2016 15:13:52 -0800*

* More work on xpbin.py (closing issues #42 and #52).


*ximpol (0.17.0) - Tue, 02 Feb 2016 15:41:14 -0800*

* Major refactoring of xpbin.py (issue #42).
* Minimum and maximum valid times added to the model components.
* Configuration file for a GRB added.


*ximpol (0.16.1) - Tue, 26 Jan 2016 18:49:19 -0800*

* Minor refactoring of the ximpol.core.fitsio module.


*ximpol (0.16.0) - Tue, 26 Jan 2016 18:40:11 -0800*

* Module ximpol.core.fitsio added (issue #49).
* ximpol.evt.event refactored to use the new ximpol.core.fitsio module.
* GTI list in the output event file (issue #24)
* ROI source table in the output event file (issue #45).
* IRF name added in the output event file header (issue #24).
* ROI information added in the output event file header (issue #48).


*ximpol (0.15.2) - Mon, 25 Jan 2016 18:04:33 -0800*

* Minor refactoring of bin/xpimgview.py


*ximpol (0.15.1) - Mon, 25 Jan 2016 16:37:52 -0800*

* astropy.wcs used in ximpol/srcmodel/img.py, and aplpy still used for
  plotting (issue #41).
* Documentation for ximpol/srcmodel/img.py added.


*ximpol (0.15.0) - Mon, 25 Jan 2016 15:57:27 -0800*

* srcmodel config files renamed.
* Point source in the Crab complex sample file dimmer.
* Added option to xpimgview.py to save the image to file.
* Horrible hack in the azimuthal fit to prevent the visibility from going
  negative (issue #34, significantly more work needed).
* Some refactoring and more documentation.
* Radius removed from the xROIModel class, and ROI model for the Crab
  nebula now correctly centered on the right coordinates.


*ximpol (0.14.0) - Fri, 22 Jan 2016 20:54:23 -0800*

* xpobbsim.py generating an output file name based on the source model
  (if not specified).
* Added CMAP mode to xpbin.py


*ximpol (0.13.0) - Fri, 22 Jan 2016 13:58:51 -0800*

* Implemented the infrastructure for multiple source in ROI


*ximpol (0.12.1) - Fri, 22 Jan 2016 06:44:01 -0800*

* Bug fix in srcmodel/source.py.


*ximpol (0.12.0) - Thu, 21 Jan 2016 16:35:14 -0800*

* First implementation of extended sources.


*ximpol (0.11.1) - Wed, 20 Jan 2016 16:57:24 -0800*

* Minor addition to the doc.


*ximpol (0.11.0) - Wed, 20 Jan 2016 15:43:39 -0800*

* load_irf moved from bin/xpobssim.py to irf/__init__.py, so that it can be
  reused.
* Unit test for IRF plotting added (issue #30).
* Some documentation for the IRFs added.


*ximpol (0.10.1) - Tue, 19 Jan 2016 16:41:33 -0800*

* More documentation and unit tests.


*ximpol (0.10.0) - Tue, 19 Jan 2016 14:45:50 -0800*

* Added math support in the sphinx config file.
* Major refactoring of the classes related to the modulation factor (issue #28).
* More unit tests added.
* More documentation added.


*ximpol (0.9.1) - Sat, 16 Jan 2016 07:17:52 -0800*

* All unit tests fixed (issue #26).


*ximpol (0.9.0) - Fri, 15 Jan 2016 16:34:58 -0800*

* IRFs extended ("by hand") down below 1 keV (need to do it properly, see
  issue #19).
* A couple of subtle bug fixes in the energy dispersion (see issues #21 and
  #22).
* First version that allows to recover the spectral parameters in XSPEC.


*ximpol (0.8.0) - Fri, 15 Jan 2016 11:53:01 -0800*

* Obsolete files removed, and some name refactoring.
* xpbin.py created.
* All figures from unit tests moved to doc/figures.
* More unit tests.
* Event times in xpobbsim sorted.
* Spectral analysis in xspec added.


*ximpol (0.7.0) - Thu, 14 Jan 2016 15:15:44 -0800*

* Modulation factor generator returning angles in degrees.
* Unit test for the modulation factor classes added.
* Source configuration moved out of xpobsim.py
* Folder srcmodel/config created.
* Added optimization step for the x grid in
  xInterpolatedBivariateSplineLinear.build_vppf() (issue #18).


*ximpol (0.6.3) - Wed, 13 Jan 2016 16:16:38 -0800*

* .travis.yml file tweaked to add display support for matplotlib.


*ximpol (0.6.2) - Wed, 13 Jan 2016 16:11:55 -0800*

* One more unit test added.


*ximpol (0.6.1) - Wed, 13 Jan 2016 15:38:20 -0800*

* Parameter tweak in the xEnergyDispersionMatric class.
* Added unit test for the xCountSpectrum class, with inline images.
* One unit test relaxed.


*ximpol (0.6.0) - Wed, 13 Jan 2016 12:13:06 -0800*

* Number of XIPE energy channels changed from 1024 to 256 and IRFs
  regenerated.
* Removed all the hard-coded values for the number of energy channels
  (issue #13).
* xEnergyDispersionMatrix now inheriting from xUnivariateAuxGenerator (i.e.,
  it has facilities to throw random numbers.)
* Down-sampling mechanism implemented for the xEnergyDispersionMatrix class
  on the energy axis to streamline performance.


*ximpol (0.5.0) - Tue, 12 Jan 2016 15:24:17 -0800*

* A couple of bug fixes in the irf.mrf module.
* Major xpobbsim refactoring using all the new classes.


*ximpol (0.4.2) - Mon, 11 Jan 2016 07:08:21 -0800*

* Minor refactoring.


*ximpol (0.4.1) - Sun, 10 Jan 2016 08:01:03 -0800*

* Grid optimization for the spline definition implemented (issue #15).
* Small application for visualizing an event file (xpevtview.py) created,
  and plotting stuff moved out of xpobbsim.


*ximpol (0.4.0) - Sat, 09 Jan 2016 10:17:52 -0800*

* New module ximpol.core.rand created (issue #16).
* Major rework and speed up of the provisional observation simulator (event
  loop removed).
* New event list classe in.
* Some cleanup.


*ximpol (0.3.1) - Thu, 07 Jan 2016 16:36:04 -0800*

* Added PSF classes, with facility to draw random numbers.


*ximpol (0.3.0) - Thu, 07 Jan 2016 13:53:07 -0800*

* Added make_ppf to the spline base class.
* Some improvement in the plotting facility for the energy dispersion.
* Added unit tests for the irf classes.
* Removed the xmin and xmax arguments from the constructor of all the spline
  classes, since the integral() method does not understand extrapolations and
  having spurious values outside the array ranges was causing troubles.
  (Note the splines can still be extrapolates in the evaluation.)
* Added facilities for normalization, cdf and ppf in the univariate spline
  base class.
* xmerge() method of the base univariate spline class removed in favor of
  numpy.union1d()


*ximpol (0.2.1) - Thu, 07 Jan 2016 06:57:12 -0800*

* First full implementation of the energy dispersion.


*ximpol (0.2.0) - Wed, 06 Jan 2016 15:56:38 -0800*

* Refactoring of the core.spline module, and plotting functionalities added.
* Unit tests for the utils.os module added.
* Initial import of the utils.matplotlib configuration module.
* Added xEffectiveArea class to irf.arf.
* Added xModulation factor class to mrf.arf.
* bin/xpirfview.py application created (issue #7).


*ximpol (0.1.2) - Tue, 05 Jan 2016 08:34:30 -0800*

* Minor changes.


*ximpol (0.1.1) - Tue, 05 Jan 2016 07:05:43 -0800*

* Minor refactoring of the irf specifications, with the OGIP part now included
  in ximpol.irf.base
* Some documentation added to the irf classes.


*ximpol (0.1.0) - Mon, 04 Jan 2016 16:15:30 -0800*

* setup.py file added (issue #11).
* release folder renamed as tools.
* ximpol.__logging__ module moved to ximpol.utils.logging (issue #8).
  Note we use the trailing undescore to avoid name conflicts with the
  correponding module from the standard library.)
* ximpol.__utils__ module splitted into ximpol.utils.os and
  ximpol.utils.system (issue #8).
* Code to create the instrument response functions moved to detector.xipe.
* New spline code used when generating the response functions and old
  xFunction1d classes removed (issue #3).
* fileio folder removed.
* Using the astropy facilities to generate the fits headers (issue #4).


*ximpol (0.0.16) - Sun, 03 Jan 2016 14:31:56 -0800*

* ximpol is now linked to Travis CI, and the build output is shown and linked
  from the main github page.


*ximpol (0.0.15) - Sat, 02 Jan 2016 07:19:39 -0800*

* xChrono class moved to utils.profile. Documentation and unit tests in place.


*ximpol (0.0.14) - Sat, 02 Jan 2016 06:59:19 -0800*

* Minor formatting fix.


*ximpol (0.0.13) - Sat, 02 Jan 2016 06:56:54 -0800*

* Added a makefile for the unit tests, and some more documentation about it.


*ximpol (0.0.12) - Fri, 01 Jan 2016 07:51:56 -0800*

* Some more edits and additions to the documentation.
* Module core.xInterpolatedUnivariateSpline moved to core.spline.
* __package__.py removed, and content moved to ximol.__init__.py, with all
  imports changed accordingly (issue #10).
* Code to be executed in __main__ moved from test() to main() in all modules
  (since the test code will be in the form of unit tests).


*ximpol (0.0.11) - Thu, 31 Dec 2015 17:19:37 -0800*

* Started migrating the documentation from the github wiki to the rst sphinx
  files, and added more stuff.


*ximpol (0.0.10) - Wed, 30 Dec 2015 07:53:08 -0800*

* Bug fix in the release script (hopefully).


*ximpol (0.0.9) - Wed, 30 Dec 2015 07:48:26 -0800*

* Major folder restructuring to make the layout compatible with
  `Read the Docs <https://readthedocs.org/>`_.
* Documentation effort started (issue #1).
* Suite of unit tests started (issue #4).
* These release notes moved to a .rst file (issue #12).
* utils.xFunction1d being replaced by core.xInterpolatedUnivariateSpline


*ximpol (0.0.8) - Mon, 28 Dec 2015 06:29:54 -0800*

* Added script to generate the rmf file. Still not working perfectly.
* Some folder refactoring.


*ximpol (0.0.7) - Fri, 11 Dec 2015 13:33:49 -0800*

* Removed the srcmodel/yaml folder and all the associated parser classes.


*ximpol (0.0.6) - Fri, 11 Dec 2015 06:39:21 -0800*

* Many minor changes.
* First stab at a parser for the source model.
* FITS images of some sources added, along with a small visualization script.
* Added a script that generates the header for the mrf file.
* Added a script to generate the .mrf file based on the ascii table provided.


*ximpol (0.0.5) - Tue, 08 Dec 2015 11:41:24 -0800*

* Small fix in the .arf XIPE file.


*ximpol (0.0.4) - Tue, 08 Dec 2015 11:33:40 -0800*

* Added a first stab at the effective area table definition.
* Added ascii data files for the XIPE IRFs (as in the proposal).
* Script to generate the .arf file for XIPE based on the ascii table.
* Added a general-purpose one-dimensional function class.


*ximpol (0.0.3) - Fri, 04 Dec 2015 12:11:49 -0800*

* Changed thge release note because I was cheating...


*ximpol (0.0.2) - Fri, 04 Dec 2015 12:05:42 -0800*

* Folder structure created


*ximpol (0.0.1) - Fri, 04 Dec 2015 06:39:19 -0800*

* Initial setup of the github repository

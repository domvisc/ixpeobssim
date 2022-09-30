from __future__ import print_function, division
from tkinter import N

import os

from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_

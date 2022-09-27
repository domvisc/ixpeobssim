from xspec import *
import os
from matplotlib import pyplot as plt

from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII


os.chdir(IXPEOBSSIM_CONFIG_ASCII)

s_A = Spectrum('nu80802313002_A_opt')
s_B = Spectrum('nu80802313002_B_opt')

#m_1 = Model('tbabs*diskbb',setPars=(5.5,1.48,149.))
m_2 = Model('tbabs*(diskbb+powerlaw)',setPars=(6.,1.45,166.,2.6,0.26))
#m_3 = Model('tbabs*bbody',setPars=(7,1))

Fit.perform()
Fit.show()
Plot.device = "/xs"
Plot.xAxis = "channel"
Plot.add = True #shows sum models
Plot('data','resid')

chans = Plot.x()
rates = Plot.y()
folded = Plot.model()

plt.plot(chans, rates, '.', chans, folded)
plt.xlabel('channels, model = tbabs*(diskbb+plaw)')
plt.ylabel('counts/cm^2/sec/chan')
plt.show()


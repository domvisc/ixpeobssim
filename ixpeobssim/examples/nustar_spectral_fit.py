from xspec import *
import os
from matplotlib import pyplot as plt

from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII


os.chdir(IXPEOBSSIM_CONFIG_ASCII)

s_A = Spectrum('nu80802313002_A_opt')
s_B = Spectrum('nu80802313002_B_opt')

#m_1 = Model('tbabs*diskbb',setPars=(5.5,1.48,149.))
m_2 = Model('tbabs*(diskbb+powerlaw)',setPars=(5.8,1.45,164.,2.6,0.25))
#m_3 = Model('tbabs*bbody',setPars=(7,1))
#m_4 = Model('tbabs*kerrbb',setPars =(0.,0.,60,4.,1.,15.,1.5,-1.,-1.,1.))

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


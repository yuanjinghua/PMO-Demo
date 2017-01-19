import aplpy as ap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
import astropy.coordinates as coords

mpl.rc("font", family="serif", size=12)
mpl.rc("axes", linewidth = 1.0)
mpl.rc("lines", linewidth = 1.0)
mpl.rc("xtick.major", pad = 5, width = 1)
mpl.rc("ytick.major", pad = 5, width = 1)
mpl.rc("xtick.minor", width = 1)
mpl.rc("ytick.minor", width = 1)


fitsDir = '../fitsDir/'
sourList = ascii.read('../sourList.txt')

epsDir = '../epsDir/'
fitsDir = '../fitsDir/'


fig2 = plt.figure(figsize = (5,5))
fb1  = ap.FITSFigure(fitsDir+'G171.3+02.5_01_13co_m1.fits', 
	figure = fig2)
fb1.show_colorscale(vmin=-20.5, vmax = -18.1, cmap = 'jet')
fb1.show_colorbar(axis_label_text = 'km s$^{-1}$')
fb1.tick_labels.set_xformat("hhmmss")
fb1.tick_labels.set_yformat("ddmm")
fb1.show_beam()
levs = np.linspace(0.3,0.9,7)*22 #*np.max(m0_13co.hdu.data)
fb1.show_contour(fitsDir+'G171.3+02.5_01_13co_m0.fits', 
	levels = levs, colors = 'red')
fb1.add_label(0, 1.03, 'G171.3+02.5_01', relative = True, 
              horizontalalignment='left', color = "black")
fb1.add_label(1, 1.03, '$^{13}$CO moment 1', relative = True, 
              horizontalalignment='right', color = "black")

fig2.savefig(epsDir+'G171.3+02.5_01_m1.pdf', bbox_inches='tight',
		papertype='a2')

fig2 = plt.figure(figsize = (5,5))
fb1  = ap.FITSFigure(fitsDir+'G171.3+02.5_01_13co_m2.fits', 
	figure = fig2)
fb1.show_colorscale(vmin= 0.1, vmax = 1.3, cmap = 'jet')
fb1.show_colorbar(axis_label_text = 'km s$^{-1}$')
fb1.tick_labels.set_xformat("hhmmss")
fb1.tick_labels.set_yformat("ddmm")
fb1.show_beam()
levs = np.linspace(0.3,0.9,7)*22 #*np.max(m0_13co.hdu.data)
fb1.show_contour(fitsDir+'G171.3+02.5_01_13co_m0.fits', 
	levels = levs, colors = 'red')
fb1.add_label(0, 1.03, 'G171.3+02.5_01', relative = True, 
              horizontalalignment='left', color = "black")
fb1.add_label(1, 1.03, '$^{13}$CO moment 2', relative = True, 
              horizontalalignment='right', color = "black")

fig2.savefig(epsDir+'G171.3+02.5_01_m2.pdf', bbox_inches='tight',
		papertype='a2')
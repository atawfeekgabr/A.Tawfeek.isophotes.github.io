import numpy as np           # to define our table
#import math                  # for computing the maximun ellipticity and fbar
#import decimal               # for computing the maximun ellipticity and fbar
#import fitsio                # for fits image
from astropy.io import fits  # isteade of fitsio
import aplpy
from astropy.io import fits
from astropy import units as u
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.ticker import MaxNLocator
from pylab import *
from photutils.isophote import EllipseGeometry, Ellipse         # main package for isophotes and ellipticity
from photutils import EllipticalAperture
from photutils.isophote import build_ellipse_model
import sys

scale=0.2

pi=22/7

#pa = pa * np.pi/180

hdu=fits.open('A3128_G183_B.fits')
datafit = hdu[0].data
fig1 = aplpy.FITSFigure(hdu[0])

g = EllipseGeometry(x0=98, y0=98, sma=10, eps=0.65, pa=12* np.pi/180,fix_center=True)
ell = Ellipse(datafit, geometry=g)
isolist = ell.fit_image()
print("ell_len", len(isolist))
sma= isolist.sma
print("sma_len", len(sma))

plt.imshow(datafit, origin='lower')
fig1.show_grayscale(vmin=None, vmid=None, vmax=None, pmin=0.25, pmax=99.75)
fig1.add_scalebar(15.*u.arcsec)
fig1.scalebar.set_corner('bottom right')
fig1.scalebar.set_label('15"') # is not applied if coming before corner
fig1.scalebar.set_color('white')
fig1.scalebar.set_linewidth('3')
fig1.scalebar.set_font_size('20')
fig1.set_yaxis_coord_type('scalar')
fig1.set_xaxis_coord_type('scalar')
#fig1.tick_labels.hide()
#fig1.axis_labels.hide()
#fig1.set_title('A3128_G64_B')
#fig1.set_theme('publication')

smas = np.linspace(1, 35, 10)
for sma in smas:
    iso = isolist.get_closest(sma)
    x, y, = iso.sampled_coordinates()
    plt.plot(x, y, color='red', alpha=0.5)
    #plt.imshow(hdu[0].data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
#plt.MaxNLocator(nbins=20, steps=[1,2,3,4,5,6,7,8,9,10])
plt.contour(hdu[0].data, levels=20, colors='yellow', alpha=0.5)
plt.text(5,180,'G183_B',color='white',fontsize=20)


#aper = EllipticalAperture((g.x0, g.y0), g.sma, g.sma * (1 - g.eps),g.pa)
#aper.plot(color='white')
#plt.imshow(datafit, origin='lower')
#data_image = build_ellipse_model(datafit.shape, isolist)
#plt.imshow(data_image, origin='lower')
#plt.show()

fig1.save('A3128_G183_B.pdf')

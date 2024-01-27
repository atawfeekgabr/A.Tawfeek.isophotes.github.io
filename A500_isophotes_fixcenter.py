import numpy as np           # to define our table
import math                  # for computing the maximun ellipticity and fbar
import decimal               # for computing the maximun ellipticity and fbar
#import fitsio                # for fits image
from astropy.io import fits  # isteade of fitsio
import bz2
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pylab import *
from photutils.isophote import EllipseGeometry, Ellipse         # main package for isophotes and ellipticity 
from photutils import EllipticalAperture
import sys

scale=0.2

pi=22/7                      # used in fbar equation     #fbar= 2/pi [arctan(1-emax)^(-1/2) - arctan(1-emax)^(1/2)]




# First define your table for geometric parameters

d= "A500_centroid_pixel.dat"
data = np.loadtxt(d,  dtype={'names':('Index','RA','DEC','WINGS','Cluster','SB','x0', 'y0', 'sma', 'eps', 'pa'),'formats': ('i','O','O','O','O', 'f','f', 'f', 'f', 'f', 'f')}, skiprows=1)

# then define what each column has

Index=data['Index']
RA= data['RA']
Dec=data['DEC']
WINGS=data['WINGS']
Cluster=data['Cluster']
SB=data['SB']
x0 = data['x0']
y0 = data['y0']
sma = data['sma']
eps = data['eps']
pa = data['pa']*np.pi/180   # we have to convert the PA into radians


# to print the output values in a new table
f_out = open('A500_output_fixcenter.dat', 'w')
print ("#Index RA Dec WINGS Cluster SB x0  y0  sma  eps  pa edelta1 edelta2 fbar1 fbar2 Bar", file=f_out)  

# Define your image (cluster)


hdu=fits.open('A500_V.fits')
datafit = hdu[0].data
#img = ma.masked_array(img, img==0)
for row in data:
    [Index,RA, Dec, WINGS, Cluster, SB, x0, y0, sma, eps, pa] = row
    pa = pa * np.pi/180
    print ("Iter %s..." % Index)  #print the object name before its output isophote table

    g = EllipseGeometry(x0, y0, sma, eps, pa, fix_center=True)

    ell = Ellipse(datafit, geometry=g)

    isolist = ell.fit_image()
    #print(isolist.to_table())
    print("E_all_len", len(isolist.eps)) 
    E_all=isolist.eps 
    sma_all= isolist.sma
    print("sma_all_len", len(sma_all)) 
    sma_max=max(sma_all)
    print('max_sma=', sma_max)
    
    
    E1_li= []       # corresponding values of E in the first third of sma
    sma1_li=[]      #for the first third of the sma
    smabar1_li= []  # sma region after the first Emax
    Ebar1_li= []    # E region after Enmax1 to identify the drop (Emin)
    sma2_li=[]      #for the second third of the sma
    E2_li= []       # corresponding values of E in the second third of sma
    smabar2_li= []  # sma region after the second Emax
    Ebar2_li= []    # E region after Emax2
    
    
    for i in range(len(sma_all)):
        if 3 < sma_all[i] <=(max(sma_all)*0.33):  ## 0.33 is 1/3 the length of the sma
            sma_1=sma_all[i] 
            #print('sma1', sma_1)
            sma1_li.append(sma_1)
            #print('sma1_li', sma1_li)
            E1= E_all[i]
            #print('E1:', E1)
            E1_li.append(E1)
            #print('E1_li', E1_li)
            E1_array=np.array(E1_li)
            
                   
        if (max(sma_all)*0.33) <=  sma_all[i] <= (max(sma_all)*0.67):   ## 0.67 is 2/3 the length of the sma
            sma_2=sma_all[i]
            #print('sma2', sma_2)
            sma2_li.append(sma_2)
            #print('sma2_li', sma2_li)
            E2= E_all[i]
            #print('E2:', E2)
            E2_li.append(E2)
            #print('E2_li', E2_li)
            E2_array=np.array(E2_li)
            
            
            
    Emax1= np.max(E1_array)
    print('Emax1:', Emax1)
    sma1_x = sma1_li[np.argmax(E1_array)]
    print('sma1_x',sma1_x)
    #print('sma1_li_len', len(sma1_li))
    #print('E1_len', len(E1_li))
    Emax2= np.max(E2_array)
    print('Emax2:', Emax2)
    sma2_x = sma2_li[np.argmax(E2_array)]  # to estimate the value of sma corresponding to Emax 
    #sma_x=np.where(E2_array==Emax)        # to estimate the number of array (e.g [13]) that has sma corresponding to Emax 
    #print('sma2_x',sma2_x)
    
    mag=22.5 -2.5*log(isolist.intens/3631) + 2.5*log(0.04)
    
    for i in range(len(sma_all)):
        if sma1_x <= sma_all[i]  <= (sma_max*0.33):
            sma1_bar=sma_all[i]
            smabar1_li.append(sma1_bar)
            #print('sma1_bar',smabar1_li)
            Ebar1=E_all[i]
            Ebar1_li.append(Ebar1)
            #print('Ebar1_li', Ebar1_li)
            #print(type (Ebar1_li))
            Emin1=min(Ebar1_li)
                
                    
        if sma2_x <= sma_all[i] <= (sma_max*0.67):
            sma2_bar=sma_all[i]
            smabar2_li.append(sma2_bar)
            #print('sma2_bar',smabar2_li)
            Ebar2=E_all[i]
            Ebar2_li.append(Ebar2)
            #print('Ebar2_li', Ebar2_li)
            Emin2=min(Ebar2_li)
                
    Edelta1=Emax1-Emin1
    Edelta2=Emax2-Emin2
    
    if Edelta1 >= 0.1:
        fbar1= (2/pi)*(np.arctan(1-(Emax1))**(-1/2) - np.arctan(1-(Emax1))**(1/2))
        flag1=round(fbar1,2)
        print('fbar1:', fbar1)
    elif Edelta1 < 0.1:
        flag1='NF'
        #print("NF")
    
    if Edelta2 >= 0.1:
        fbar2= (2/pi)*(np.arctan(1-(Emax2))**(-1/2) - np.arctan(1-(Emax2))**(1/2))
        flag2=round(fbar2,2)
        print('fbar2:', fbar2)
    elif Edelta2 < 0.1:
        flag2='NF'
        #print("NF")
        
    if flag1 == flag2 == 'NF':
        Bar='F' 
    else: 
        Bar='T'
                               
        
    print(Index, RA, Dec, WINGS, Cluster, SB, x0,  y0,  sma,  eps, round(pa/np.pi*180.,1), round(Edelta1,2),round(Edelta2,2), flag1, flag2, Bar, file=f_out)
    
       
    
    plt.figure(1)
    plt.plot(isolist.sma*scale, mag, 'o' )
    plt.gca().invert_yaxis()
    plt.xlabel('Semimajor axis ($arcsec$)')
    plt.ylabel('SB ($mag/arcsec^2$)')
    savefig('A500_SB_G%s.png' %Index)
    plt.close()
    
    plt.figure(2)
    fig, axs = plt.subplots(2)
    axs[0].errorbar(isolist.sma, isolist.eps, yerr=isolist.ellip_err, fmt='o', markersize=4)
    axs[0].set( xlabel='Semimajor axis ($pixels$)', ylabel='Ellipticity')
    axs[1].errorbar(isolist.sma, isolist.pa/np.pi*180., yerr=isolist.pa_err/np.pi* 80., fmt='o', markersize=4)
    axs[1].set(xlabel='Semimajor axis ($pixels$)', ylabel='PA (deg)')
    savefig('A500_eps_G%s.png' %Index)
    plt.close
    
    


                

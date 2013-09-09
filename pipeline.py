#from reduce_obs import *
from getztest import *
from multigauss import *
from lineprops import *
from coplot import *
from scipy import integrate

import numpy as np
import matplotlib.pyplot as plt
import os

def coprops(source,start=1,wtf=2,nsims=100000,g=1.,D=32.,H0=[70,2],wm=[.27,.01],wv=[.73,.01],zk=0,s1100=0.,s14=0.,e1100=0.,e14=0.):
    '''Pipeline for reducing RSR data. 
1. Search for all observations with name starting with source
2. Cross correlation method with rejection based on absence of a bright line
3. Estimate amplitude, frequency, and fwhm based on the redshift from 2
4. Cross correlation method with convolution of gaussian and fwhm from 3 for
   a better estimate for the redshift.
5. Estimate initial conditions from the new best fit redshift
6. Find the best fit and errors in the fit from the initial conditions in 5
7. From best fit parameters find z, fwhm, Sco, D, and Lco
Inputs
source: SourceName
start:  If 0, reduce observations, else if a txt file for source exists, skip
        step 1
wtf:    Which weighting function to use for x-corr method. 1 (default) = temp  
        late, 0= equal weighting, 2 = bright lines
nsims:  Number of sims for line fitting technique
H0:     Hubble parameter value/error
wm:     omega matter value/error
wv:     omega vacuum value/error
zk:     if the redshift of an object is known, redshift input
Outputs
xcorz:   column1: redshifts, column2:cross correlation values
freqspec:column1: frequency, column2:spectrum, column3:convolved spectrum
         column4: best fit spectrum
lp:      array1: best fit redshift, array2: distance, array3: fwhm,
         array4: SCO, array5: Lco, array6: initial conditions for fit
    '''
    ospath = './%s.txt' % source
#    if os.path.exists(ospath): #checks if the data file exists, if it does and the reduction is done properly, head into the analysis steps
#        start = start
#    else:
#        start = 0
    
#    if start == 0:
#        reduce_obs(source) #gives the option to analyze the data if the file doesn't exist or the reduction needs to be redone

    fname = '%s.txt' % source
    z,bf,zs,x,y,xcor,z1,z21 = getz(filename = np.loadtxt(fname),wtarray=wtf,s1100=s1100,s14=s14,e1100=e1100,e14=e14) #cross correlation to determine if there is a best redshift
    dz = zs[np.argsort(zs)]
    fs = np.transpose([x,y])
    ics,nlines = initcons(filename=fs,zbf=z,g=g,D=D)
    if np.max(fs[:,1])<10*np.std(fs[:,1]):
        fs1 = cnvlvspec(filename=fs,sig=.05) #convolve the spectrum with a gaussian of FWHM ~1000 km/s if there are no strong spectral lines 
    else:
        fs1 = fs
    zbf2,bf,zs,x,y,xcor,z1,z21 = getz(filename=fs1,wtarray=wtf,s1100=s1100,s14=s14,e1100=e1100,e14=e14) #redo cross correlation if there are no strong spectral lines
    ics,nlines = initcons(filename=fs,zbf=zbf2,g=g,D=D) #determine the initial conditions for the MCMC fit based on the peak value in the spectrum and the redshift    
    x1,y,ybf,abf,ebf = multigauss(filename=fs,ics=ics,nsims=nsims,nlines=nlines,g=g,D=D) #multi-line MCMC fit 
    zbf =  np.zeros((2,nlines)) #preallocate the physical and dynamic properties of the galaxy
    fwhm = np.zeros((2,nlines))
    Sco = np.zeros((2,nlines))
    Dl = np.zeros((2,nlines))
    Lco = np.zeros((2,nlines))
    if zk != 0:
        zbf2 = zk

    for j in range(nlines):
        zbf[:,j],fwhm[:,j],Sco[:,j],Dl[:,j],Lco[:,j] = lineprops(fname,z=zbf2,abf=abf[j,:],ebf=ebf[j,:],H0=H0,wm=wm,wv=wv) #line properties based on the best fit parameters and the spectrum
        Sco = Sco/1000. #assigning the proper physical units to the derived line properties
        Lco10 = Lco/10**10
        zprint= 'z = %.4f +- %.4f' % (zbf[0,j], zbf[1,j])
        if z21>=.5 or zk != 0:
            dprint= 'D_L = %.0f +- %.0f Mpc' % (Dl[0,j], Dl[1,j])
        else:
            dprint = 'Redshift Ambiguous: Distance Unknown'    
        fprint= 'fwhm = %.0f +- %.0f km/s' % (fwhm[0,j], fwhm[1,j])
        sprint= 'S_CO = %.1f +- %.1f Jy km/s' % (Sco[0,j], Sco[1,j])
        if z21>=1. or zk != 0:
            lprint= 'L_CO = %.2f +- %.1f 10^10  K km/s pc' % (Lco10[0,j], Lco10[1,j])
        else:
            lprint = 'Redshift Ambiguous: Luminosity Unknown'

        print 'Observations of ', source #print all these values
        print zprint
        print dprint
        print fprint
        print sprint
        print lprint

    xbf = np.arange(73,111,.001)
    coplot(fs[:,0],y,xbf,ybf,zbf[0,0],source=source) #create publication quality figures of the observed source
    xcorz = np.transpose([dz,xcor])
    lp = [zbf,fwhm,Sco,D,Lco,ics]
    return lp #returns the best fit parameters

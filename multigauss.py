import scipy as sp
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
#import idlsave

def initcons(filename,zbf,g=1.,D=30.):
    '''
    Blindly determine initial guesses for the MCMC fit for the source
    based on the most likely redshift obtained in getz.py
    '''
    x1 = filename[:,0]
    Ta = filename[:,1]*1000
    lam = 3e10/(x1*10**9)
    theta = lam/30e2*206264.8
    #y1 = 1.08*(x1/115.)**2*theta**2/100.*Ta/.5
    y1 = 4654.*Ta*g*D**-2
    freqwe = np.genfromtxt('lines.catalog')
    lfreq1 = freqwe[:,0]
    wt1 = freqwe[:,3]
    linename = np.loadtxt('lines.catalog',dtype=str)
    name11 = linename[:,1]
    name21 = linename[:,2]
    wt = []
    name1 = []
    name2 = []
    lfreq = []
    for j in range(np.size(wt1)):
        if wt1[j] == 1 and lfreq1[j]/(1+zbf)>73 and lfreq1[j]/(1+zbf)<111:
            wt = np.append(wt,wt1[j])
            name1 = np.append(name1,name11[j])
            name2 = np.append(name2,name21[j])
            lfreq = np.append(lfreq,lfreq1[j])

    x = []
    y = []
    for j in range(np.size(x1)):
        if y1[j] < np.inf:
            x = np.append(x,x1[j])
            y = np.append(y,y1[j])

    sval = np.zeros(np.size(y)/2)
    sfre = np.zeros(np.size(x)/2)
    for j in range(np.size(sval)):
        sval[j] = np.average(y[2*j:2*(j+1)])
        sfre[j] = np.average(x[2*j:2*(j+1)])
    nlines = np.size(wt)
    ics = np.zeros((3,nlines))
    for j in range(nlines):
        ics[0,j] = lfreq[j]/(1+zbf)
        ics[1,j] = np.max(y)
        mline = np.argsort(sval)[np.size(sval)-1]
        dn = 0
        for k in range(1,5):
            if sval[mline+k]<.5*np.max(sval) and sval[mline+k-1]>.5*np.max(sval) and dn == 0:
                ics[2,j] = (sfre[2]-sfre[1])*k/(np.sqrt(2*np.log(2)))
                dn = 1
    return ics,nlines
    
    

def multigauss(filename,ics,nlines,nsims=100000,g=1.,D=30.):
    '''
    Gaussian MCMC fit of spectral lines observed in each source
    '''
    x1 = filename[:,0]
    Ta = filename[:,1]*1000
    lam = 3e10/(x1*10**9)
    theta = lam/30e2*206264.8
    y1 = 4654.*Ta*g*D**-2
    errnum = np.ceil(nsims*.683)
    x = []
    y = []
    for j in range(np.size(x1)):
        if y1[j] < np.inf:
            x = np.append(x,x1[j])
            y = np.append(y,y1[j])

    if np.shape(ics)[0] != 3 or np.shape(ics)[1] != nlines:
        print 'Initial conditions should have shape (3,nlines)'
    else:
        a0v = np.zeros((nlines,nsims))
        a1v = np.zeros((nlines,nsims))
        a2v = np.zeros((nlines,nsims))
        abf = np.zeros((nlines,3))
        ebf = np.zeros((nlines,3))
        v = np.zeros((nlines,1))
        ev = np.zeros((nlines,1))
        step = []
        chis = np.zeros((nlines,nsims))
        ybf = 0.
        for j in range(nlines):
            print j
            xm1 = (x-ics[0,j])**2.
            xm = np.argsort(xm1)[0]
            totd = np.ceil(2*np.sqrt(2*np.log(ics[1,j]/np.std(y)/.01))*ics[2,j]/.03125)
            #print xm
            if xm<totd:
                x2 = x[:(xm+totd)]
            elif xm>(np.size(x)-totd):
                x2 = x[(xm-totd):]
            else:
                x2 = x[(xm-totd):(xm+totd)]
            a = ics[:,j]
            f = lambda a,x2: a[1]*np.exp(-(x2-a[0])**2/(2*a[2]**2))
            for k in range(nsims):
                yn = y[(xm-totd):(xm+totd)]
                ya = f(a,x2)
                step = np.append(step,k)
                b = [0,0,0]
                b[0] = a[0]+np.random.randn(1)*.01
                b[1] = a[1]+a[1]*np.random.randn(1)*.015*np.std(y)
                b[2] = np.abs(a[2]+a[2]*np.random.randn(1)*.003)
                if k<1000.: #burn in period for MCMC fit..roughly find the best fit parameters in a better way than the initcons function
                    b[0] = a[0]+np.random.randn(1)*.01
                    b[1] = a[1]+a[1]*np.random.randn(1)*.1*np.std(y)
                    b[2] = np.abs(a[2]+a[2]*np.random.randn(1)*.05)
                yb = f(b,x2)
                piba = 1
                for l in range(np.size(x2)):
                    piba = piba*np.exp(ya[l]**2-2.*yn[l]*ya[l]+2.*yn[l]*yb[l]-yb[l]**2)
                alpha = np.min((1.,piba))
                #print piba
                u = np.random.rand(1)
                if u <= alpha:
                    a = b
                else:
                    a = a
                chisq = np.sum((yn-f(a,x2))**2.)
                chis[j,k] = chisq
                a0v[j,k] = a[0]
                a1v[j,k] = a[1]
                a2v[j,k] = a[2]
            chisbf = np.argsort(chis[j,:])[0]
            abf[j,:] = [a0v[j,chisbf],a1v[j,chisbf],a2v[j,chisbf]]
            mm = (abf[j,0]-a0v[j,1000:])**2.
            mm1 = np.argsort(mm)[errnum]
            ma = (abf[j,1]-a1v[j,1000:])**2.
            ma1 = np.argsort(ma)[errnum]
            ms = (abf[j,2]-a2v[j,1000:])**2.
            ms1 = np.argsort(ms)[errnum]
            ebf[j,:] = [np.abs(a0v[j,1000+mm1]-abf[j,0]),np.abs(a1v[j,1000+ma1]-abf[j,1]),np.abs(a2v[j,1000+ms1]-abf[j,2])]
            ebf[j,0] = np.sqrt(ebf[j,0]**2)
            ebf[j,1] = np.sqrt(ebf[j,1]**2)
            ebf[j,2] = np.sqrt(ebf[j,2]**2)
            xbf = np.arange(73,111,.001)
            ybf = ybf+abf[j,1]*np.exp(-(xbf-abf[j,0])**2./(2*abf[j,2]**2))
    return x, y, ybf, abf, ebf

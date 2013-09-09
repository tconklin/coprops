import scipy as sp
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import idlsave

def photoz(s1100,e1100=0.,s14=0.,e14=0.,ntry=50000):
    '''
    Determine the photometric redshift of a galaxy given the
    measured 1.4 cm and 1100 micron flux and uncertainty
    '''
    z = np.arange(0,10,.05)
    ngal = 44
    if s14 == 0:
        ratioin = -1
        ratiosig = -1
    else:
        ratioin = s1100/s14
        ratiosig = (e1100/s1100**2+e14/s14**2)**.5
    a = idlsave.read('fluxratio1100.sav')
    dat = a.get('data')
    zs = a.get('redshift')
    averatio = np.zeros(200)
    sigma = np.zeros(200)
    array = np.random.randn(ntry)
    array1 = np.random.randn(ntry)
    if s14 <= 0.:
        ydarts = (s1100+array*e1100)/(np.abs(array1*e14))
    else:
        ydarts = array*ratiosig+ratioin
    xdarts = np.zeros(ntry)
    for i in range(ntry):
        jrangal = np.floor(ngal*np.random.rand(1))[0]
        testtrack = dat[:,jrangal]
        yval = ydarts[i]
        xdarts[i] = np.interp(yval,testtrack,z)
    return xdarts,ydarts

def rmnan(filename):
    '''
    Remove nan values from the spectrum (this occurs because the band
    edges of a spectrometer typically have much worse response). Nan
    values cannot be used to in determining statistics
    '''
    sfreq1 = filename[:,0]
    svalue1 = filename[:,1]
    ndata = np.size(sfreq1)
    #print ndata, ' lines read from the input spectrum file'
    #print ''
    svalue=[]
    sfreq=[]
# I don't want any "nan"
    for j in range(ndata):
        if svalue1[j]<np.inf:
            svalue = np.append(svalue,svalue1[j])
            sfreq = np.append(sfreq,sfreq1[j])
    return sfreq,svalue


def redshiftgen(freqspec,sn=10,z1=0,z2=10):
    '''
    generates a random spectrum with arbitrary signal to noise ratio
    between two redshifts. This was done including the actual response
    of the instrument, however I just use gaussian noise in this 
    github repository
    '''
    z = np.random.rand(1)*(z2-z1)+z1
    #print z
    sfreq,svalue = rmnan(freqspec)
    sortval = np.argsort(sfreq)
    svalue = np.ones(np.size(sortval))
    svalue0 = svalue
    sfreq = sfreq[sortval]
    freqwe = np.genfromtxt('lines.catalog')
    lfreq = freqwe[:,0]
    wt = freqwe[:,3]
    linename = np.loadtxt('lines.catalog',dtype=str)
    name1 = linename[:,1]
    name2 = linename[:,2]
    nlines = np.size(lfreq)
    wtarray = np.zeros((4,np.size(sfreq)))
    lfreqz = lfreq/(1+z)
    for i in range(nlines):
        if lfreqz[i]>np.min(sfreq) and lfreqz[i]<np.max(sfreq):
            ichannum = np.argsort(np.abs(lfreqz[i]-sfreq))[0]
            wtarray[:,ichannum]=wtarray[:,ichannum]+wt[i]
            if ichannum>3 and ichannum<(np.size(sfreq)-4):
                wtarray[1,ichannum-1]=wtarray[1,ichannum-1]+.8*wt[i]
                wtarray[1,ichannum+1]=wtarray[1,ichannum+1]+.8*wt[i]
                wtarray[2,ichannum-1]=wtarray[2,ichannum-1]+.8*wt[i]
                wtarray[2,ichannum+1]=wtarray[2,ichannum+1]+.8*wt[i]
                wtarray[2,ichannum-2]=wtarray[2,ichannum-2]+.7*wt[i]
                wtarray[2,ichannum+2]=wtarray[2,ichannum+2]+.7*wt[i]
                wtarray[3,ichannum-1]=wtarray[3,ichannum-1]+.8*wt[i]
                wtarray[3,ichannum+1]=wtarray[3,ichannum+1]+.8*wt[i]
                wtarray[3,ichannum-2]=wtarray[3,ichannum-2]+.7*wt[i]
                wtarray[3,ichannum+2]=wtarray[3,ichannum+2]+.7*wt[i]
                wtarray[3,ichannum-3]=wtarray[3,ichannum-3]+.5*wt[i]
                wtarray[3,ichannum+3]=wtarray[3,ichannum+3]+.5*wt[i]
    wtarray = wtarray[3,:]*sn
    wts = wtarray*svalue0+np.random.randn(np.size(svalue0))
    fs = np.zeros((np.size(sfreq),2))
    fs[:,0] = sfreq
    fs[:,1] = wts
    return fs,z

def weights(filename,z1=0,z2=10,dz=.001,wtf=0):
    '''
    Determine the weights of the cross correlation matrix. A variety
    of different weighting schemes such as equal weighting, bright
    spectral lines, (anticipated) strength weighting are possible.
    Pre-allocating the space for these directories reduces the analysis
    time by about a factor of 100.
    '''
    x = filename[:,0]
    y = filename[:,1]
    sfreq1 = filename[:,0]
    svalue1 = filename[:,1]
    ndata = np.size(sfreq1)
    print ndata, ' lines read from the input spectrum file'
    print ''
    svalue=[]
    sfreq=[]
    sfreq,svalue=rmnan(filename)
    sortval = np.argsort(sfreq)
    svalue = svalue[sortval]
    sfreq = sfreq[sortval]
    sval = np.zeros(np.size(svalue)/1)
    sfre = np.zeros(np.size(svalue)/1)
    for j in range(np.size(sval)):
        sval[j] = np.average(svalue[1*j:1*(j+1)])
        sfre[j] = np.average(sfreq[1*j:1*(j+1)])

    svalue = sval
    sfreq = sfre 
    z = np.arange(z1,z2,dz)
    nz = np.size(z)
    wtarray  = np.zeros((nz,4,np.size(sfreq)))
    freqwe = np.genfromtxt('lines.catalog')
    lfreq = freqwe[:,0]
    wt = freqwe[:,3]
    linename = np.loadtxt('lines.catalog',dtype=str)
    name1 = linename[:,1]
    name2 = linename[:,2]
    nlines = np.size(lfreq)
    print nlines, ' lines read from the catalog'
    print ''
    if wtf==0: #equal weighting all lines
        for j in range(nz):
            lfreqz = lfreq/(1.+z[j])
            for i in range(nlines):
                zfreq = lfreqz[i]
                if zfreq>np.min(sfreq) and zfreq<np.max(sfreq):
                    ichannum = np.argsort(np.abs(zfreq-sfreq))[0]
                    wtarray[j,:,ichannum]=wtarray[j,:,ichannum]+1
                    if ichannum>3 and ichannum<(np.size(sfreq)-4):
                        #print ichannum
                        wtarray[j,1,ichannum-1]=wtarray[j,1,ichannum-1]+1
                        wtarray[j,1,ichannum+1]=wtarray[j,1,ichannum+1]+1
                        wtarray[j,2,ichannum-1]=wtarray[j,2,ichannum-1]+1
                        wtarray[j,2,ichannum+1]=wtarray[j,2,ichannum+1]+1
                        wtarray[j,2,ichannum-2]=wtarray[j,2,ichannum-2]+1
                        wtarray[j,2,ichannum+2]=wtarray[j,2,ichannum+2]+1
                        wtarray[j,3,ichannum-1]=wtarray[j,3,ichannum-1]+1
                        wtarray[j,3,ichannum+1]=wtarray[j,3,ichannum+1]+1
                        wtarray[j,3,ichannum-2]=wtarray[j,3,ichannum-2]+1
                        wtarray[j,3,ichannum+2]=wtarray[j,3,ichannum+2]+1
                        wtarray[j,3,ichannum-3]=wtarray[j,3,ichannum-3]+1
                        wtarray[j,3,ichannum+3]=wtarray[j,3,ichannum+3]+1
        np.save('wtf0.npy',wtarray)
    elif wtf==1: # template spectra weighting
        for j in range(nz):
            #print j
            lfreqz = lfreq/(1.+z[j])
            for i in range(nlines):
                #print i
                zfreq = lfreqz[i]
                if zfreq>np.min(sfreq) and zfreq<np.max(sfreq):
                    ichannum = np.argsort(np.abs(zfreq-sfreq))[0]
                    #print wtarray[j,:,ichannum]
                    wtarray[j,:,ichannum]=wtarray[j,:,ichannum]+wt[i]
                    if ichannum>3 and ichannum<(np.size(sfreq)-4):
                        #print ichannum
                        wtarray[j,1,ichannum-1]=wtarray[j,1,ichannum-1]+wt[i]
                        wtarray[j,1,ichannum+1]=wtarray[j,1,ichannum+1]+wt[i]
                        wtarray[j,2,ichannum-1]=wtarray[j,2,ichannum-1]+wt[i]
                        wtarray[j,2,ichannum+1]=wtarray[j,2,ichannum+1]+wt[i]
                        wtarray[j,2,ichannum-2]=wtarray[j,2,ichannum-2]+wt[i]
                        wtarray[j,2,ichannum+2]=wtarray[j,2,ichannum+2]+wt[i]
                        wtarray[j,3,ichannum-1]=wtarray[j,3,ichannum-1]+wt[i]
                        wtarray[j,3,ichannum+1]=wtarray[j,3,ichannum+1]+wt[i]
                        wtarray[j,3,ichannum-2]=wtarray[j,3,ichannum-2]+wt[i]
                        wtarray[j,3,ichannum+2]=wtarray[j,3,ichannum+2]+wt[i]
                        wtarray[j,3,ichannum-3]=wtarray[j,3,ichannum-3]+wt[i]
                        wtarray[j,3,ichannum+3]=wtarray[j,3,ichannum+3]+wt[i]
        np.save('wtf1.npy',wtarray)
    elif wtf==2: # strong line weighting
        for j in range(nz):
            lfreqz = lfreq/(1.+z[j])
            for i in range(nlines):
                zfreq = lfreqz[i]
                if zfreq>np.min(sfreq) and zfreq<np.max(sfreq) and wt[i]>.2:
                    ichannum = np.argsort(np.abs(zfreq-sfreq))[0]
                    wtarray[j,:,ichannum]=wtarray[j,:,ichannum]+wt[i]
                    if ichannum>3 and ichannum<(np.size(sfreq)-4):
                        #print ichannum
                        wtarray[j,1,ichannum-1]=wtarray[j,1,ichannum-1]+wt[i]
                        wtarray[j,1,ichannum+1]=wtarray[j,1,ichannum+1]+wt[i]
                        wtarray[j,2,ichannum-1]=wtarray[j,2,ichannum-1]+wt[i]
                        wtarray[j,2,ichannum+1]=wtarray[j,2,ichannum+1]+wt[i]
                        wtarray[j,2,ichannum-2]=wtarray[j,2,ichannum-2]+wt[i]
                        wtarray[j,2,ichannum+2]=wtarray[j,2,ichannum+2]+wt[i]
                        wtarray[j,3,ichannum-1]=wtarray[j,3,ichannum-1]+wt[i]
                        wtarray[j,3,ichannum+1]=wtarray[j,3,ichannum+1]+wt[i]
                        wtarray[j,3,ichannum-2]=wtarray[j,3,ichannum-2]+wt[i]
                        wtarray[j,3,ichannum+2]=wtarray[j,3,ichannum+2]+wt[i]
                        wtarray[j,3,ichannum-3]=wtarray[j,3,ichannum-3]+wt[i]
                        wtarray[j,3,ichannum+3]=wtarray[j,3,ichannum+3]+wt[i]
        np.save('wtf2.npy',wtarray)
    return wtarray

def cnvlvspec(filename,sig=.05):
    '''
    Convolve the spectrum with a gaussian with an arbitrary width.
    '''
    freq,spec = rmnan(filename)
    x1 = np.arange(-.2,.2,freq[1]-freq[0])
    y1 = 1/(sig*np.sqrt(np.pi*2))*np.exp(-.5*(x1/sig)**2)
    print np.trapz(y1,x1)
    cspec = np.convolve(spec,y1,'same')*(freq[1]-freq[0])
    cspec *= np.max(spec)/np.max(cspec)
    fs = np.transpose([freq,cspec])
    return fs

def getz(filename,z1=0.,z2=10.,dz=.001,outfile=' ',wtarray=2,sourcename=' ',s1100=0.,s14=0.,e1100=0.,e14=0.):
    '''
    Cross correlation analysis of the spectrum. Statistically determine
    the most likely redshift of the source. A paper is in preparation
    describing this technique and how it can be used.
    '''
    sfreq2,svalue2 = rmnan(filename)
    sortval = np.argsort(sfreq2)
    svalue2 = svalue2[sortval]
    sfreq2 = sfreq2[sortval]
    sval = np.zeros(np.size(svalue2)/1)
    sfre = np.zeros(np.size(svalue2)/1)
    for j in range(np.size(sval)):
        sval[j] = np.average(svalue2[1*j:1*(j+1)])
        sfre[j] = np.average(sfreq2[1*j:1*(j+1)])

    svalue = sval
    sfreq = sfre 

# remove a DC offset with poly_fit
    nfit = 2
    yval = np.polyfit(sfreq,svalue,nfit)
    yfit = yval[0]*sfreq**2+yval[1]*sfreq+yval[2]
    svalue0 = svalue - yfit
    svalsort = np.argsort(svalue0)
    svalstd = svalue0[svalsort]
    sfreqstd = sfreq[svalsort]
    stdsval = np.std(svalstd[0:(np.size(svalstd)-30)]) 

# Instrument details
    freq0 = np.min(sfreq)        # freq at channel 1
    bfreq = freq0       # beginning of the band
    efreq = np.max(sfreq)   # end of the band

# get the line freq. and  weight array from the line catalog
    freqwe = np.genfromtxt('lines.catalog')
    lfreq = freqwe[:,0]
    wt = freqwe[:,3]
    linename = np.loadtxt('lines.catalog',dtype=str)
    name1 = linename[:,1]
    name2 = linename[:,2]
    nlines = np.size(lfreq)

# loop through the redshift and calculate the cross-corr func.
    z = np.arange(z1,z2,dz)
    nz = np.size(z)
    xcor = np.zeros((4,nz))
    xcorztxt = 'xcorz%s.npy' % (wtarray)
    xcorz = np.load(xcorztxt)
    if np.size(wtarray)<=2:
        if wtarray==0:
            wtarray = np.load('wtf0.npy')
        elif wtarray==1:
            wtarray = np.load('wtf1.npy')
        elif wtarray==2:
            wtarray = np.load('wtf2.npy')

    for j in range(nz):
        #print j
        wts = np.zeros(np.shape(svalue0))
        #print np.shape(svalue0), ' ', np.shape(wtarray[j,:,:])
        wts = wtarray[j,:,:]*svalue0
        sumwts = np.sum(wtarray[j,0,:])
        #print np.max(wts)
        xcor[0,j]=np.sum(wts[0,:])
        #print xcor[0,j]
        xcor[1,j]=np.sum(wts[1,:])
        xcor[2,j]=np.sum(wts[2,:])
        xcor[3,j]=np.sum(wts[3,:])
        blines = np.argsort(wtarray[j,0,:])[(np.shape(wtarray)[2]-5):]

    xcor = xcor*xcorz/np.max(xcorz)
    if s1100>0.:
        xd,yd = photoz(s1100,e1100,s14,e14)
        histz,dhistz,listc = plt.hist(xd,normed='True',range=[0,6],bins=24)
        histz2 = np.zeros(np.size(dhistz))
        histz2[1:] = histz
        xcorh = np.interp(z,dhistz,histz2)
        xcor *= xcorh
    xcorder = np.argsort(xcor)
    zs = z[xcorder[1,:]]
# ###########################################
# Find the maximum correlation amplitude in z
#
    maxval = np.max(xcor,1)
    zbf = z[xcorder[1,nz-1]]
    lfreqzb = lfreq/(1+zbf)
    bestarray = np.zeros(np.size(sfreq))

    for i in range(nlines):    
        if lfreqzb[i]>bfreq and lfreqzb[i]<efreq:
            ichannum = np.argsort(np.abs(lfreqzb[i]-sfreq))[0]
            bestarray[ichannum]=bestarray[ichannum]+wt[i]

    bestfit = bestarray*np.max(svalue0)
    vc = v = 3.e5*(zbf*(zbf+2.)/(2.+2.*zbf+zbf**2))
    
    zbftext = 'Best Fit z = %s' % (zbf)
    xa = xcor[1,:][np.argsort(xcor[1,:])]
    alpha0 = (z2-z1)/dz-1
    alpha = alpha0
    while np.abs(zs[alpha0]-zs[alpha])<.01:
        alpha = alpha-1
    z1 = xa[alpha0]
    z21 = (xa[alpha0]-xa[alpha])/(np.std(xa))
    return zbf,bestfit,zs,sfreq,svalue,xcor[1,:],z1,z21

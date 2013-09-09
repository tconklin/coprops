import scipy as sp
import numpy as np
from scipy import interpolate
from scipy import integrate
from scipy import stats
import matplotlib.pyplot as plt
#import idlsave

def lineprops(filename,z,abf,ebf,D=32,H0=[70,2],wm=[.27,.01],wv=[.73,.01]):
    '''
    Determine the best fit physical and intrinsic conditions inside
    the observed galaxy. These are determined through scientific 
    principles and using other radiative transfer codes (not developed
    by me)
    '''
    c = 3.e5
    k = 1.38e16
    c1 = c*10**5
    freqwe = np.genfromtxt('lines.catalog')
    lfreq1 = freqwe[:,0]
    wt1 = freqwe[:,3]
    linename = np.loadtxt('lines.catalog',dtype=str)
    fs = np.loadtxt(filename)
    f = fs[:,0]
    sp = fs[:,1]*4654*1./D**2*1000
    snu = stats.nanstd(sp)
    nc = np.ceil(2*np.sqrt(2*np.log(1000))*abf[2]/.031)
    fm = np.argsort((f-abf[0])**2)[:nc]
    fm = fm[np.argsort(fm)]
    name11 = linename[:,1]
    name21 = linename[:,2]
    wt = []
    name1 = []
    name2 = []
    lfreq = []
    for j in range(np.size(wt1)):
        if wt1[j] == 1:
            wt = np.append(wt,wt1[j])
            name1 = np.append(name1,name11[j])
            name2 = np.append(name2,name21[j])
            lfreq = np.append(lfreq,lfreq1[j])

    ebf[0] = (ebf[0]**2+(.031/np.sqrt(8*np.log(2)))**2)**.5
    ebf[2] = (ebf[2]**2+(.031/np.sqrt(8*np.log(2)))**2)**.5
    rfreq = lfreq/(1+z)
    nlines = np.size(abf)/3.
    bfline = lfreq[np.argsort((abf[0]-rfreq)**2)[0]]
    zbf = [(bfline/abf[0]-1),(bfline/abf[0]**2*ebf[0])]
    #print bfline
    #v = (c*(1+z)**2-1)/(1+(1+z)**2)
    zs = np.linspace(0,zbf[0],1000)
    d = c/H0[0]*np.trapz((wm[0]*(1+zs)**3+wv[0])**-.5,zs)
    D = [d*(1+zbf[0]),0]
#    eDwm = c/(H0[0])*(1+zbf[0])*np.trapz(((wm[0]+wm[1])*(1+zs)**3+wv[0]-wm[1])**-.5,zs)-D[0]
    D[1] = D[0]*H0[1]/H0[0]
    fwhm1 = [2*np.sqrt(2*np.log(2))*abf[2]/abf[0]*c,0]
    fwhm = [fwhm1[0]-2*np.sqrt(2*np.log(2))*(fwhm1[0]-(fwhm1[0]**2-(.031*c/abf[0])**2)**.5),0]
    R = c*.031/abf[0]
    fwhm[0] += -70.*np.sqrt(8.*np.log(2.))*(np.sqrt(1+(R/70.)**2/(8.*np.log(2.)))-1.)
    fwhm[1] = ((ebf[0]*c/abf[0])**2+(ebf[2]*c/abf[0])**2)**.5
    Sco = [integrate.simps(sp[fm],f[fm])*c/abf[0],0]
    Sco[1] = np.sqrt(3*fwhm[0]*c*.031/abf[0])*snu
    Lco = [2.350*Sco[0]*(115/abf[0])**2*D[0]**2*(1+zbf[0])**-3,0] #Sco in mJy
    Lco[1] = Lco[0]*((Sco[1]/Sco[0])**2+(2*D[1]/D[0])**2+(2*ebf[0]/abf[0])**2)**.5
    return zbf,fwhm,Sco,D,Lco

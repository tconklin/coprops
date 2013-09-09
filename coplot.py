import numpy as np
import matplotlib.pyplot as plt

def coplot(x,y,xbf=0,ybf=0,zbf=0,source=''):
    '''
    Plotting utility for the Redshift Search Receiver. Create
    publication quality images of each spectral line observed. The
    final image will show a zoomed in view of the spectral line with
    the entire spectrum shown as an inset
    '''
    plt.rcParams.update({'font.size':22}) #set up the plotting window
    freqwe = np.genfromtxt('lines.catalog') #read in the spectral line catalog
    lfreq = freqwe[:,0]
    wt = freqwe[:,3]
    linename = np.loadtxt('lines.catalog',dtype=str)
    name1 = linename[:,1]
    name2 = linename[:,2]
    nlines = np.size(lfreq)
    lfreqzb = lfreq/(1+zbf) #doppler shift the spectrum to what is observed
    for i in range(nlines):    
        if lfreqzb[i]>np.min(x) and lfreqzb[i]<np.max(x) and wt[i]>.8: #plot the lines that are observed
            ichannum = np.argsort(np.abs(lfreqzb[i]-x))[0]
            name3 = '%s (%s)' % (name1[i],name2[i]) #set up the annotations for the figures and plot the spectrum
            plt.close()
            plt.figure(figsize=(25,10))
            plt.plot(x,y,linestyle='steps')
            plt.plot(xbf-.0165,ybf)
            plt.annotate(name3,xycoords='data',xy=(x[ichannum],y[ichannum]*1.1),size='14',rotation='vertical')
            plt.xlabel('Frequency (GHz)')
            plt.ylabel('Flux (mJy)')
            bestz = 'z = %.3f' % zbf
            plt.annotate(bestz,xycoords='axes fraction',xy=(.9,.9)) #plot the inset of the entire spectrum
            titled = '%s %s' % (source, name3)
            plt.title(titled)
            plt.axis([lfreqzb[i]-4,lfreqzb[i]+1,np.min(y),1.4*np.max(y)])
            plt.axes([.2,.6,.4,.2])
            plt.plot(x,y,linestyle='steps')
            plt.xlabel('Frequency (GHz)',size='14')
            plt.ylabel('Flux (mJy)')
            fname = '%s%s%s.png' % (name1[i],name2[i],source)
            plt.savefig(fname,dpi=200) #save the figure as a png file
            plt.show()

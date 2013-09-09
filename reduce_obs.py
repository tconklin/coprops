from dreampy.redshift.netcdf import RedshiftNetCDFFile
from dreampy.utils.filterscans import FilterScans
from dreampy.redshift.plots import RedshiftPlot
import os
import numpy as np #import required modules
def reduce_obs(sourceobs,thresh_sig=.1,interactive=True):
    '''
    Reduce all observations of an observed source. 
    '''
    fs = FilterScans() #set up filter scans
    fs.filter(Receiver='RedshiftReceiver',ObsPgm='Bs',Source=sourceobs) #filter scans based on IRAS
    obslist = [] # obslist
    el = []
    for j in range(10000): #how many times to run the loop
        try:
            h = fs.observations[j] #holder for the jth observation
            obslist = np.append(obslist,h.obsnum) #appends obslist
        except:
            continue
    obslist = np.unique(obslist)
    pl = RedshiftPlot()
    count = 0 #sizee of sig/temp arrays
    sig=[] #empty array for sigma values
    temp = [] #empty array for temp values
    ospath = './%s.nc' % sourceobs
    if os.path.exists(ospath):
        os.remove(ospath) #deletes old file if exists
    ncout = RedshiftNetCDFFile(ospath, 'w') #new nc4 file in write mode
    tint = 0.0 #integration time, currently 0
    #what this loop does: gets an obsnum, cycles through all chassis, if the obsnum 
    for ObsNum in obslist: #for observations in obslist
        for chassis in (0, 1, 2): #fol all chassis
            try:
                if ObsNum in (3267,3304,3306,3321,3624,3633,3673,3791,3802,3824,3833,3903,3982,4079,4092): #determined a priori, there is no way of knowing which scans are bad
                    continue
                fname = fs.get_redshift_filename(ObsNum=ObsNum, ScanNum=1,
                                                 chassis=chassis) # gets file name for current scan in obsnum, chassis
                if fname:
                    print "Process filename %s" % fname
                    nc = RedshiftNetCDFFile(os.path.join('/raw_mnc',fname))
                else:
                    continue
            except:
                continue
            if ObsNum == 3957 and chassis == 1:
                continue
            print nc.hdu.header.SourceName
            count += 1
            nc.hdu.process_scan()
            nc.hdu.baseline(order=0, subtract=False)
            nc.hdu.average_all_repeats(weight='sigma')
            temp=np.append(temp,nc.hdu.header.get('Rsfend.Temperature'))
            sig=np.append(sig,nc.hdu.sigma)
            integ = 2*int(nc.hdu.header.get('Bs.NumRepeats'))*int(nc.hdu.header.get('Bs.TSamp'))
            el = np.append(el,nc.hdu.header.get('Sky.ElReq')[0]*180./3.14159)
            tint += integ
            if interactive == True:
                pl.clf()
                pl.plot_spectra(nc)
                zz = 1
                zz = raw_input('To reject observation, type ''n'':')
                if zz != 'n':
                    nc.append_observation(ncout)
            else:
                if np.mean(sig)<thresh_sig:
                    nc.append_observation(ncout)
            nc.close()
            del nc

    print tint
    sig = np.reshape(sig,(count,6))
    ncout.close()
    objnc = '%s.nc' % sourceobs
    ncobj = RedshiftNetCDFFile(objnc)
    grp,objhdu=ncobj.average_all_hdus(threshold_sigma=thresh_sig)
    #pl = RedshiftPlot()
    if interactive==True:
        pl.plot_spectra(objhdu)
        baselinesub = raw_input('Order of baseline (type ''n'' for none):')
        if baselinesub == 'n':
            objhdu.baseline(order=0, subtract=False)
        elif baselinesub == '':
            objhdu.baseline(order=0, subtract=True)
        else:
            objhdu.baseline(order=int(baselinesub),subtract=True)
    else:
        objhdu.baseline(order=0,subtract=True)
    txtfl = '%s.txt' % sourceobs
    objhdu.make_composite_scan()
    objhdu.write_composite_scan_to_ascii(txtfl)

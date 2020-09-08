import sys

import numpy as np
from astroquery.mast import Observations, Catalogs
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib

import subprocess

import json

import tfit5 #FORTRAN based transit model
import math #needed for floor command

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib
except ImportError:  # Python 2.x
    import httplib


def transitplot(time,flux,sol,nplanetplot=1, itime=-1, ntt=0, tobs=0, omc=0, dtype=0):
    "plot the transit model"
    nplanet=int((len(sol)-8)/10) #number of planets

    #deal with input vars and translating to FORTRAN friendly.
    
    if type(itime) is int :
        if itime < 0 :
            itime=np.ones(len(time))*0.020434
        else:
            itime=np.ones(len(time))*float(itime)
    
    if type(ntt) is int :
        nttin=  np.zeros(nplanet, dtype="int32") #number of TTVs measured 
        tobsin= np.zeros(shape=(nplanet,len(time))) #time stamps of TTV measurements (days)
        omcin=  np.zeros(shape=(nplanet,len(time))) #TTV measurements (O-C) (days)
    else:
        nttin=ntt
        tobsin=tobs
        omcin=omc
 
    if type(dtype) is int :
        dtypein=np.ones(len(time), dtype="int32")*int(dtype) #contains data type, 0-photometry,1=RV data
    
    #remove other planets for plotting
    sol2=np.copy(sol)
    for i in range(1,nplanet+1):
        if i!=nplanetplot:
            nc=8+10*(i-1)
            sol2[nc+3]=0.0 #rdrs
    tmodel= np.zeros(len(time)) #contains the transit model
    tfit5.transitmodel(nplanet,sol2,time,itime,nttin,tobsin,omcin,tmodel,dtypein)
    stdev=np.std(flux-tmodel)
    #print(stdev)
    
    #make a model with only the other transits to subtract
    nc=8+10*(nplanetplot-1)
    sol2=np.copy(sol)
    sol2[nc+3]=0.0 #rdrs
    tmodel2= np.zeros(len(time)) #contains the transit model
    tfit5.transitmodel(nplanet,sol2,time,itime,nttin,tobsin,omcin,tmodel2,dtypein)
    
    epo=sol[nc+0] #time of center of transit
    per=sol[nc+1] #orbital period
    zpt=sol[7] #photometric zero-point
    
    tdur=tfit5.transitdur(sol,1)/3600.0 #transit duration in hours
    
    ph1=epo/per-math.floor(epo/per) #calculate phases
    phase=[]
    tcor=tfit5.lininterp(tobsin,omcin,nplanetplot,nttin,epo)
    #print(tcor,nttin,tobsin[1,1],omcin[1,1])
    for x in time:
        if nttin[nplanetplot-1] > 0:
            tcor=tfit5.lininterp(tobsin,omcin,nplanetplot,nttin,x)
        else:
            tcor=0.0
        t=x-tcor
        ph=(t/per-math.floor(t/per)-ph1)*per*24.0 #phase in hours offset to zero.
        phase.append(ph)
        
    phase = np.array(phase) #convert from list to array
    
    phasesort=np.copy(phase)
    fluxsort=np.copy(tmodel)
    p=np.ones(len(phase), dtype="int32") #allocate array for output. FORTRAN needs int32 
    tfit5.rqsort(phase,p)    
    i1=0
    i2=len(phase)-1
    for i in range(0,len(phase)):
        phasesort[i]=phase[p[i]-1]
        fluxsort[i]=tmodel[p[i]-1]
        if phasesort[i] < -tdur:
            i1=i
        if phasesort[i] <  tdur:
            i2=i
    
    fplot=flux-tmodel2+1.0
    
    plt.figure(figsize=(12,10)) #adjust size of figure
    matplotlib.rcParams.update({'font.size': 22}) #adjust font
    plt.scatter(phase, fplot, c="blue", s=100.0, alpha=0.35, edgecolors="none") #scatter plot
    plt.plot(phasesort, fluxsort, c="red", lw=3.0)
    plt.xlabel('Phase (hours)') #x-label
    plt.ylabel('Relative Flux') #y-label
    x1,x2,y1,y2 = plt.axis()    #get range of plot
    ymin=min(fluxsort[i1:i2])
    ymax=max(fluxsort[i1:i2])
    y1=ymin-0.1*(ymax-ymin)-2.0*stdev
    y2=ymax+0.1*(ymax-ymin)+2.0*stdev
    plt.axis((-tdur,tdur,y1,y2)) #readjust range
    plt.show()  #show the plot
        
    return;

def build_initial_sol(starpars,exocat,toi_index):

    #Build solution
    sol=[]
    serr=[]

    #mean-stellar density
    sol.append(starpars.rhostar)
    serr.append(-1.0) #fitted parameter

    #limb-darkening

    sol.append(0.0) #u1,u2
    sol.append(0.0)
    serr.append(0.0)
    serr.append(0.0)
    #Use q1, q2 formalism
    sol.append(starpars.q1)
    sol.append(starpars.q2)
    serr.append(-1.0)
    serr.append(-1.0)

    #dil,vof
    sol.append(0.0) #Dilution 
    serr.append(0.0)
    sol.append(0.0) #radial Velocity offset
    serr.append(0.0)

    #zpt
    sol.append(0.0) #photometric offset
    serr.append(-1)

    #t0
    sol.append(exocat.t0[toi_index])
    serr.append(-1.0)
    #per
    sol.append(exocat.per[toi_index])
    serr.append(-1.0)
    #b
    sol.append(0.4) #Impact parameter is missing from NEXSCI table.
    serr.append(-1.0)
    #rdr 
    sol.append(np.sqrt(exocat.tdep[toi_index]/1.0e6)) #also a guess, since no r/R* we guess from Tdepth
    serr.append(-1.0)

    #ecw,esw #These work well with e~<0.99
    sol.append(0.0) #sqrt(e)cos(w)
    serr.append(0.0)
    sol.append(0.0) #sqrt(e)sin(w)
    serr.append(0.0)

    #KR,TE,EL,AL
    sol.append(0.0) #radial velocity amplitude
    serr.append(0.0) 
    sol.append(0.0) #Secondary Eclipse Depth
    serr.append(0.0)
    sol.append(0.0) #Ellipsodial variations (scales by a^3)
    serr.append(0.0)
    sol.append(0.0) #Albedo amplitude (Lambertian Sphere Approx.)
    serr.append(0.0)
    
    return sol,serr


class starpars_class:
    def __init__(self):
        self.teff=[]  #initialize arrays
        self.teff_err=[]
        self.logg=[]
        self.logg_err=[]
        self.feh=[]
        self.feh_err=[]
        rhostar=[]
        rhostar_err=[]
        rhostar_ep=[]
        rhostar_em=[]
        q1=[]
        q2=[]
        q1_err=[]
        q2_err=[]

def get_limb_q1q2(teff,logg,feh,bindir='./bin/'):
    """Quick hack that calls a command line executable to get limb-darkening coefficients.
    This will be replaced by a proper module. 
    """
    
    
    output=subprocess.check_output([bindir+"claretquad_tess",\
                                    str(teff),str(logg),str(feh),str(1)]).decode('ascii')
    #print(output)
    for row in output.split('\n'):
        nline=row.split()
        if len(nline) == 6:
            #titles.append(nline[0])
            u1=np.float(nline[0])
            u2=np.float(nline[1])
            
            q1=np.power(u1+u2,2)
            q1=max(min(1.0,q1),0.0) # 0 <=q1<=1
            q2=u1/(2*(u1+u2))
            q2=max(min(1.0,q2),0.0) # 0 <=q2<=1
            
    return q1,q2

def limbprior(starpars,nsamp=100):
    """Get Limb-darkening priors
    
    Currently assumes sysmetric errors.  
    
    starpars - starpars_class()
    """
    
    q1,q2=get_limb_q1q2(starpars.teff,starpars.logg,starpars.feh)
    
    q1samp=[]
    q2samp=[]
    for i in range(nsamp):
        teff1=starpars.teff+np.random.normal()*starpars.teff_err
        logg1=starpars.logg+np.random.normal()*starpars.logg_err
        feh1=starpars.feh+np.random.normal()*starpars.feh_err
        q1_1,q2_1=get_limb_q1q2(teff1,logg1,feh1)
        q1samp.append(q1_1)
        q2samp.append(q2_1)
        #print(q1_1,q2_1,teff1,logg1,feh1)
        
    q1samp=np.array(q1samp)
    q2samp=np.array(q2samp)
    q1sig=np.std(q1samp)
    q2sig=np.std(q2samp)
        
    
    return q1,q2,q1sig,q2sig

def mastQuery(request):

    server='mast.stsci.edu'

    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
        "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)

    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content

def ticAdvancedSearch(id):
    request = {"service":"Mast.Catalogs.Filtered.Tic",
                "format":"json",
                "params":{
                "columns":"*",
                "filters":[
                    {"paramName":"id",
                        "values":[{"min":id,"max":id}]}]
                        #"values":[{261136679}]}]
                }}

    headers,outString = mastQuery(request)
    outData = json.loads(outString)

    return outData



def get_tess_data(u_ticid,max_flag=16,out_dir='./download'):
    """Given a TIC-ID, return time,flux,ferr,itime
    u_ticid : (int) TIC ID

    returns lc_time,flux,ferr,int_time
    """
    tic_str='TIC' + str(u_ticid)
    #out_dir='/data/rowe/TESS/download/'

    # Search MAST for TIC ID
    obs_table=Observations.query_object(tic_str,radius=".002 deg")

    # Identify TESS timeseries data sets (i.e. ignore FFIs)
    oti=(obs_table["obs_collection"] == "TESS") & \
            (obs_table["dataproduct_type"] == "timeseries")
    if oti.any() == True:
        data_products=Observations.get_product_list(obs_table[oti])
        dpi=[j for j, s in enumerate(data_products["productFilename"]) if "lc.fits" in s]
        manifest=Observations.download_products(data_products[dpi],download_dir=out_dir)
    else:
        manifest=[]

    lc_time=[]
    flux=[]
    ferr=[]
    int_time=[]
    for j in range(0,len(manifest)):
        fits_fname=str(manifest["Local Path"][j])
        print(fits_fname)
        hdu=fits.open(fits_fname)
        tmp_bjd=hdu[1].data['TIME']
        tmp_flux=hdu[1].data['PDCSAP_FLUX']
        tmp_ferr=hdu[1].data['PDCSAP_FLUX_ERR']
        tmp_int_time=hdu[1].header['INT_TIME'] + np.zeros(len(tmp_bjd))
        tmp_flag=hdu[1].data['QUALITY']

        ii=(tmp_flag <= max_flag) & (~np.isnan(tmp_flux))
        tmp_bjd=tmp_bjd[ii]
        tmp_flux=tmp_flux[ii]
        tmp_ferr=tmp_ferr[ii]
        tmp_int_time=tmp_int_time[ii]
        # Shift flux measurements
        median_flux=np.median(tmp_flux)
        tmp_flux=tmp_flux / median_flux
        tmp_ferr=tmp_ferr / median_flux
        # Append to output columns
        lc_time=np.append(lc_time,tmp_bjd)
        flux=np.append(flux,tmp_flux)
        ferr=np.append(ferr,tmp_ferr)
        int_time=np.append(int_time,tmp_int_time)

        hdu.close()

    # Sort by time
    si=np.argsort(lc_time)
    lc_time=np.asarray(lc_time)[si]
    flux=np.asarray(flux)[si]
    ferr=np.asarray(ferr)[si]
    int_time=np.asarray(int_time)[si]

    phot=phot_class()
    phot.time=np.copy(lc_time)
    phot.flux=np.copy(flux)
    phot.ferr=np.copy(ferr)
    phot.itime=np.copy(int_time)    

    phot.itime=phot.itime/(60*24) #convert minutes to days

    return phot

class phot_class:
    def __init__(self):
        self.time=[]  #initialize arrays
        self.flux=[]
        self.ferr=[]
        self.itime=[]

class exocat_class:
    def __init__(self):
        self.ticid=[]
        self.toiid=[]
        self.toiid_str=[]
        self.ra=[]
        self.dec=[]
        self.tmag=[]
        self.t0=[]
        self.t0err=[]
        self.per=[]
        self.pererr=[]
        self.tdur=[]
        self.tdurerr=[]
        self.tdep=[]
        self.tdeperr=[]

# Read 'toi_file' (from NASA EA -- new table)
def readtoicsv(toi_file):
    exocat=exocat_class() #set up class
    f=open(toi_file)

    icount=0
    for line in f:
        line = line.strip()
        row = line.split(',') #break into columns
        if row[0][0]!='#':
            #skip comments
            icount+=1
            if icount>1:
                #skip header
                exocat.ticid.append(int(row[3]))
                exocat.toiid.append(float(row[1]))
                exocat.toiid_str.append(row[1])
                exocat.ra.append(float(row[8]))
                exocat.dec.append(float(row[12]))
                exocat.tmag.append(float(row[60]))

                if row[25]=='':
                    exocat.t0.append(-1.0)
                else:
                    try:
                        exocat.t0.append(float(row[25]) - 2.457E6) # Planet Transit Midpoint Value [BJD]
                    except:
                        print(row[25])

                if row[26]=='': exocat.t0err.append(-1.0)
                else: exocat.t0err.append(float(row[26])) # Planet Transit Midpoint Upper Unc [BJD]

                if row[30]=='': exocat.per.append(-1.0)
                else: exocat.per.append(float(row[30])) # Planet Orbital Period Value [days]

                if row[31]=='': exocat.pererr.append(-1.0)
                else: exocat.pererr.append(float(row[31])) # Planet Orbital Period Upper Unc [days]

                if row[35]=='': exocat.tdur.append(-1.0)
                else: exocat.tdur.append(float(row[35])) # Planet Transit Duration Value [hours]

                if row[36]=='': exocat.tdurerr.append(-1.0)
                else: exocat.tdurerr.append(float(row[36])) # Planet Transit Duration Upper Unc [hours]

                if row[40]=='': exocat.tdep.append(-1.0)
                else: exocat.tdep.append(float(row[40])) # Planet Transit Depth Value [ppm]

                if row[41]=='': exocat.tdeperr.append(-1.0)
                else: exocat.tdeperr.append(float(row[41])) # Planet Transit Depth Upper Unc [ppm]
    f.close()

    return exocat

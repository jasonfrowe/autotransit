import sys

import numpy as np
from astroquery.mast import Observations, Catalogs
from astropy.io import fits

import json

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

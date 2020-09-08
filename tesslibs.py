import os
import numpy as np
import matplotlib.pyplot as plt
from astroquery.mast import Observations, Catalogs
from astropy.io import fits
import sys
import json
import pickle

import transitfit5 as tf

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


def get_tess_data(u_ticid,max_flag=16):
    """Given a TIC-ID, return time,flux,ferr,itime
    u_ticid : (int) TIC ID

    returns lc_time,flux,ferr,int_time
    """
    tic_str='TIC' + str(u_ticid)
    out_dir='/data/rowe/TESS/download/'

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

    phot=tf.phot_class()
    phot.time=np.copy(lc_time)
    phot.flux=np.copy(flux)
    phot.ferr=np.copy(ferr)
    phot.itime=np.copy(int_time)    

    phot.itime=phot.itime/(60*24) #convert minutes to days

    return phot

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

def get_tess_ffi(u_ticid,tmag,ffi_dir):

    phot=tf.phot_class()

    #check to see what sectors have data for target.
    tic_str='TIC' + str(u_ticid)
    obs_table=Observations.query_object(tic_str,radius=".0001 deg")
    #print(obs_table)
    sectors = np.array(obs_table[obs_table["obs_collection"]=='TESS']["sequence_number"])

    for s in sectors:
        lcsfile=ffi_dir+'tesslcs_sector_'+str(s)+'/'+'tesslcs_tmag_'+str(int(tmag))+'_'+str(int(tmag)+1)+'/'+\
        'tesslc_'+str(u_ticid)+'.pkl'
        #print(lcsfile)
        if os.path.isfile(lcsfile):
            infile = open(lcsfile,'rb')
            data=pickle.load(infile)
            #print(len(data))
            infile.close()
            if len(phot.time)==0:
                medflux=np.median(data[6])
                phot.time=np.copy(data[4]) #time
                phot.flux=np.copy(data[6]/medflux) #corrected flux
                phot.ferr=np.copy(data[8]/medflux) #flux error
                phot.itime=np.ones(len(data[4]))*0.5/24 #integration time in days
            else:
                medflux=np.median(data[6])
                phot.time=np.concatenate((phot.time,data[4]))
                phot.flux=np.concatenate((phot.flux,data[6]/medflux))
                phot.ferr=np.concatenate((phot.ferr,data[8]/medflux))
                phot.itime=np.concatenate((phot.itime,np.ones(len(data[4]))*0.5/24))

        else:
            print("Warning. File not found: ",lcsfile)

    return phot

import sys

import numpy as np
from astroquery.mast import Observations, Catalogs
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter

import subprocess

from scipy import stats #For Kernel Density Estimation

import json

import tfit5 #FORTRAN based transit model
import math #needed for floor command

#Fortran based fancy detrender. It's a bit slow, because it does automagic gap detection.
#However it's about 100x faster than its python equivilent.
import detrend5

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

def histplots(chain,sol,serr,burnin,nbin,label,colour):

    npar_sys=8 #number of system-wide parameters (e.g., limb-darkening)
    npar_pln=10 #planet specific parameters (e.g., r/R*)

    nbodies=int((len(sol)-npar_sys)/npar_pln)
    npl=npar_sys+nbodies*npar_pln

    matplotlib.rcParams.update({'font.size': 10}) #adjust font
    plt.figure(1, figsize=(20, 8*nbodies))
    nullfmt = NullFormatter()

    n1=int(nbodies+1) #added 1
    n2=int(npar_pln)

    jpar=-1
    for j in range(npar_sys):
        if np.abs(serr[j]) > 1.0e-30:
            jpar=jpar+1
            npanel=int(jpar+1)
            axScatter=plt.subplot(n1, n2, npanel)

            minx=np.min(chain[burnin:,jpar]) #range of parameter
            maxx=np.max(chain[burnin:,jpar])

            x_eval = np.linspace(minx, maxx, num=100) #make a uniform sample across the parameter range

            kde1 = stats.gaussian_kde(chain[burnin:,jpar],0.3) #Kernel Density Estimate

            #plot the histogram
            plt.hist(chain[burnin:,jpar],nbin,histtype='stepfilled', density=True, facecolor=colour[jpar],\
                     alpha=0.6)

            #overlay the KDE
            plt.plot(x_eval, kde1(x_eval), 'k-', lw=3)

            plt.xlabel(label[jpar])
            #plt.ylabel('Probability Density')
            axScatter.yaxis.set_major_formatter(nullfmt)

    #Dscale parameter
    npanel=npanel+1
    axScatter=plt.subplot(n1, n2, npanel)
    k=chain.shape[1]-1
    minx=np.min(chain[burnin:,k]) #range of parameter
    maxx=np.max(chain[burnin:,k])
    x_eval = np.linspace(minx, maxx, num=100) #make a uniform sample across the parameter range
    kde1 = stats.gaussian_kde(chain[burnin:,k],0.3) #Kernel Density Estimate
    #plot the histogram
    plt.hist(chain[burnin:,k],nbin,histtype='stepfilled', density=True, facecolor=colour[k], alpha=0.6)
    #overlay the KDE
    plt.plot(x_eval, kde1(x_eval), 'k-', lw=3)
    plt.xlabel(label[k])
    #plt.ylabel('Probability Density')
    axScatter.yaxis.set_major_formatter(nullfmt)

    for i in range(nbodies):
        for j in range(npar_sys):
            n=npar_sys+i*npar_pln+j
            #print(n,j+7)
            if np.abs(serr[n]) > 1.0e-30:
#                print(i*n1+j,j)
                jpar=jpar+1
                npanel=int((i+1)*n2+j+1)
                axScatter=plt.subplot(n1, n2, npanel)

                minx=np.min(chain[burnin:,jpar]) #range of parameter
                maxx=np.max(chain[burnin:,jpar])

                x_eval = np.linspace(minx, maxx, num=100) #make a uniform sample across the parameter range

                kde1 = stats.gaussian_kde(chain[burnin:,jpar],0.3) #Kernel Density Estimate

                #plot the histogram
                plt.hist(chain[burnin:,jpar],nbin,histtype='stepfilled', density=True, facecolor=colour[jpar],\
                     alpha=0.6)

                #overlay the KDE
                plt.plot(x_eval, kde1(x_eval), 'k-', lw=3)

                plt.xlabel(label[jpar])
                #plt.ylabel('Probability Density')
                axScatter.yaxis.set_major_formatter(nullfmt)

    plt.show()

def plotchains(chain,label,colour,burnin,psize=0.1):
    nullfmt = NullFormatter()
    matplotlib.rcParams.update({'font.size': 10}) #adjust font
    plt.figure(figsize=(12,37)) #adjust size of figure

    x=np.arange(burnin+1,len(chain)+1,1)
    npar=len(chain[0,:])
    for i in range(0,npar):
        axScatter=plt.subplot(npar, 1, i+1)
        plt.scatter(x,chain[burnin:,i],c=colour[i],s=psize,alpha=0.1)  #plot parameter a
        plt.ylabel(label[i])                   #y-label

        x1,x2,y1,y2 = plt.axis()
        y1=np.min(chain[burnin:,i])
        y2=np.max(chain[burnin:,i])
        plt.axis((x1,x2,y1,y2))


        if i < npar-1:
            axScatter.xaxis.set_major_formatter(nullfmt)

    plt.xlabel('Iteration')           #x-label

    plt.show()

def gelmanrubin(chain,burninfrac,npt):
    "Estimating PSRF"
    M=chain.shape[1]      #number of chains
    burnin=int(chain.shape[2]*burninfrac)
    N=chain.shape[2]-burnin #assuming all chains have the same size.
    npars=chain.shape[0] #number of parameters
    pmean=np.zeros(shape=(M,npars)) #allocate array to hold mean calculations 
    pvar=np.zeros(shape=(M,npars))  #allocate array to hold variance calculations

    
    for i in range(0,M):
        for j in range(0,npars):
            pmean[i,j]=np.mean(chain[j,i,burnin:]) #Generate means for each parameter in each chain
            pvar[i,j]=np.var(chain[j,i,burnin:])   #Generate variance for each parameter in each chain
    
    posteriormean=np.zeros(npars) #allocate array for posterior means
    for j in range(0,npars):
        posteriormean[j]=np.mean(pmean[:,j]) #calculate posterior mean for each parameter
        
    #Calculate between chains variance
    B=np.zeros(npars)
    for j in range(0,npars):
        for i in range(0,M):
            B[j]+=np.power((pmean[i,j]-posteriormean[j]),2)
    B=B*N/(M-1.0)    
    
    #Calculate within chain variance
    W=np.zeros(npars)
    for j in range(0,npars):
        for i in range(0,M):
            W[j]+=pvar[i,j]
    W=W/M 
    
    
    #Calculate the pooled variance
    V=(N-1)*W/N + (M+1)*B/(M*N)
    
    dof=npt-1 #degrees of freedom 
    Rc=np.sqrt((dof+3.0)/(dof+1.0)*V/W) #PSRF from Brooks and Gelman (1997)
    
    #Calculate Ru
    #qa=0.95
    #ru=np.sqrt((dof+3.0)/(dof+1.0)*((N-1.0)/N*W+(M+1.0)/M*qa))
    
    return Rc;

def autocorr_func_1d(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf

def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i

def autocorr_new(y, c=5.0):
    f = np.zeros(y.shape[1])
    for yy in y:
        f += autocorr_func_1d(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]

# Automated windowing procedure following Sokal (1989)
def auto_window(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1



def modekdestimate(chain,burnin):
    'Estimate Mode with KDE and return KDE'
    #range of data
    minx=np.min(chain[burnin:])
    maxx=np.max(chain[burnin:])
    x_eval = np.linspace(minx, maxx, num=1000)
    kde1 = stats.gaussian_kde(chain[burnin:])#,0.3)
    modeval=[]
    modekde=0
    for x in x_eval:
        if kde1(x) > modekde:
            modekde=kde1(x)
            modeval=x
    return modeval,x_eval,kde1 ;

def intperc(x,x_eval,kde1,perc=0.6827):
    'find error bounds'
    idx = (np.abs(x_eval-x)).argmin()
    kdea=np.array(kde1(x_eval))

    n=len(x_eval)

    #print(x,idx)

    i1=1
    i2=1
    intval=0.0

    j1=np.copy(idx)
    j2=np.copy(idx)
    j1old=np.copy(j1)
    j2old=np.copy(j2)
    while intval < perc:
        j1test=np.max((0,idx-i1-1))
        j2test=np.min((n-1,idx+i2+1))
        if kdea[j1test] > kdea[j2test]:
            if j1test>0:
                j1=np.copy(j1test)
                i1=i1+1
            else:
                j1=np.copy(j1test)
                j2=np.copy(j2test)
                i2=i2+1
            #print('case1')
        else:
            if j2test<n-1:
                j2=np.copy(j2test)
                i2=i2+1
            else:
                j2=np.copy(j2test)
                j1=np.copy(j1test)
                i1=i1+1
            #print('case2')

        intval=np.trapz(kdea[j1:j2],x_eval[j1:j2])
        #print(j1,j2,intval,kdea[j1test],kdea[j2test])

        #make sure we can break from loop
        if (j1 == 0) and (j2 == n-1):  #break we reach boundaries of array
            #print('break1')
            intval=1.0
        if (j1 == j1old) and (j2 == j2old): #break if stuck in loop.
            #print('break2')
            intval=1.0

        #Update old values to check we are making progress.
        j1old=np.copy(j1)
        j2old=np.copy(j2)

    #print(x_eval[j1],x_eval[j2])
    return x_eval[j1],x_eval[j2];


def checkperT0(samples,burninfrac,nthin,sol,serr):
    
    #get indices of Period and T0 from chain
    nparsol=len(sol)
    j=-1
    iT0=-1
    iPer=-1
    for i in range(nparsol):
        if np.abs(serr[i])>1.0e-10:
            j+=1
            if i==8:
                iT0=j
            if i==9:
                iPer=j
                
    #print('ii',iT0,iPer)
    
    nburnin=len(samples)*burninfrac
    burnin=int(nburnin/nthin) #correct burnin using nthin.

    chain=np.array(samples)
    chain=chain[::nthin,:] #Thin out chain.
    #print('Thin size: ',len(chain))
    #print('Burnin: ',burnin)

    if burnin > 0:
        burnin_bak=np.copy(burnin)
    else:
        burnin=np.copy(burnin_bak)
    sigcut=4
    niter=3
    for k in range(niter):
        npars=chain.shape[1]
        test=np.copy(chain[burnin:,:])
        for i in range(npars):
            nch=test.shape[0]
            #print(nch)
            mean=np.mean(test[:,i])
            std=np.std(test[:,i])
            #print(mean,std)
            test2=[]
            for j in range(nch):
                #print(test[j,i], np.abs(test[j,i]-mean),std*sigcut)
                #input()
                if np.abs(test[j,i]-mean) < sigcut*std:
                    test2.append(test[j,:])
            test=np.array(test2)
        nch=test.shape[0]
        #print("nchains:",nch)
        chain=np.copy(test)
        burnin=0
    
    
    if iT0>-1:
        mode,x_eval,kde1=modekdestimate(chain[:,iT0],burnin)
        perc1 = intperc(mode,x_eval,kde1)
        t0_ep=np.abs(perc1[1]-mode)
        t0_em=np.abs(mode-perc1[0])
    else:
        t0_ep=0.0
        t0_em=0.0
    
    if iPer>-1:
        mode,x_eval,kde1=modekdestimate(chain[:,iPer],burnin)
        perc1 = intperc(mode,x_eval,kde1)
        per_ep=np.abs(perc1[1]-mode)
        per_em=np.abs(mode-perc1[0])
    else:
        per_ep=0.0
        per_em=0.0
        
    return t0_ep,t0_em,per_ep,per_em

def readsol (filename):
    "read in transitmodel solution"
    nplanetmax=9 #maximum number of planets that an n0.dat file can handle
    nplanet=0 #count number of planets found in the solution
    solin=np.zeros(nplanetmax*10+8) #allocate array to hold parameters. init to zero.
    serrin=np.zeros(nplanetmax*10+8)
    f = open(filename, 'r')
    for line in f:
        line = line.strip() #get rid of the \n at the end of the line
        columns = line.split() #break into columns
        if columns[0][0:3]=='RHO':
            solin[0]=columns[1]
            serrin[0]=columns[3]
        elif columns[0][0:3]=='NL1':
            solin[1]=columns[1]
            serrin[1]=columns[3]
        elif columns[0][0:3]=='NL2':
            solin[2]=columns[1]
            serrin[2]=columns[3]
        elif columns[0][0:3]=='NL3':
            solin[3]=columns[1]
            serrin[3]=columns[3]
        elif columns[0][0:3]=='NL4':
            solin[4]=columns[1]
            serrin[4]=columns[3]
        elif columns[0][0:3]=='DIL':
            solin[5]=columns[1]
            serrin[5]=columns[3]
        elif columns[0][0:3]=='VOF':
            solin[6]=columns[1]
            serrin[6]=columns[3]
        elif columns[0][0:3]=='ZPT':
            solin[7]=columns[1]
            serrin[7]=columns[3]
        elif columns[0][0:2]=='EP':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+0)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='PE':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+1)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='BB':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+2)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='RD':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+3)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='EC':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+4)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='ES':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+5)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='KR':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+6)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='TE':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+7)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='EL':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+8)
            solin[j]=columns[1]
            serrin[j]=columns[3]
        elif columns[0][0:2]=='AL':
            np1=float(columns[0][2])
            if np1>nplanet:
                nplanet=np1
            j=int(10*(np1-1)+8+9)
            solin[j]=columns[1]
            serrin[j]=columns[3]
    f.close()
    #print(nplanet)
    sol=solin[0:int(nplanet*10+8)]
    serr=serrin[0:int(nplanet*10+8)]
    return sol, serr;

def detrend(phot,detrend_win,nfitp=4,nc=20,itime=0.02083,tflag=0):
    """Wrapper for Fortran based detrender.
    
    phot = (phot.class) 
        phot.time = time stamps (np.array) [days]
        phot.flux = relative flux (np.array)
    detrend_win = detrend window (np.float) [days]
    
    nfitp = polynomial order (int) 
    nc = cadence gap detection (int).  nc*itime gives the window size for gap detection
    itime = integration time (np.float) [days]
    
    """
    npt=len(phot.time)    
    if type(tflag) is int :
        tflag=np.zeros(npt,dtype=int)

    ts=np.zeros(npt,dtype=int)
    boxbin=np.copy(detrend_win)
    x=np.zeros(npt) #work arrays
    y=np.zeros(npt)
    z=np.zeros(npt)
    ngap=0
    gaps=np.zeros(npt)
    offset=np.zeros(npt)
    work=np.zeros(npt)
    
    phot_d=phot_class()
    phot_d.time=np.copy(phot.time)
    phot_d.flux=np.copy(phot.flux)
    phot_d.ferr=np.copy(phot.ferr)
    phot_d.itime=np.copy(phot.itime)
    
    detrend5.polyfilter(phot_d.time,phot_d.flux,phot_d.ferr,ts,tflag,boxbin,x,y,z,ngap,gaps,\
                        offset,nfitp,work,nc,itime)
    
    return phot_d

def write_n1dat(toi_id,sol,serr,modeltype,filedir):

    npl=int(100*(0.001+toi_id-np.floor(toi_id))) #planet number
    nfile = 'tables/nsample.dat'

    npars=0
    titles=[]
    output=subprocess.check_output(["cat",nfile]).decode('ascii')
    for row in output.split('\n'):
        nline=row.split()
        if len(nline) == 5:
            npars+=1
            titles.append(nline[0])
    
    #create n1.dat file
    nfilenew='n_'
    if modeltype[0]==1:
        nfilenew+='ld_'
    if modeltype[1]==1:
        nfilenew+='rh_'
    if modeltype[2]==1:
        nfilenew+='ec_'
    if modeltype[3]==1:
        nfilenew+='gpm32_'
    nfilenew+=str(npl)+'.dat' 
    
    
    f=open(filedir+nfilenew,'w')
    for i in range(npars):
        f.write('%s %.10e %.10e %.10e %.10e\n' %(titles[i],sol[i],0.0,serr[i],0.0))
    f.close()


def setup_walkers(sol,serr,phot,rsuf,rmin,rmax,nwalkers,modeltype):

    npar_sys=8 #number of system-wide parameters (e.g., limb-darkening)
    npar_pln=10 #planet specific parameters (e.g., r/R*)
    eps=np.finfo(float).eps #A small number
    favg=np.median(phot.ferr) #average error 

    nbodies=int((len(sol)-npar_sys)/npar_sys)
    npl=npar_sys+nbodies*npar_pln
    x=[]
    xerr=[]
    for i in range(npl): #npl only includes transit model parameters.
        if np.abs(serr[i]) > 1.0e-30:
            x.append(sol[i])

            if serr[i]<0:
                if i<npar_sys:
                    #print(i,rsuf[i],rmin[i],rmax[i])
                    xerr.append(rsuf[i])
                if i>=npar_sys:
                    j=i-int((i-npar_sys+1)/npar_pln)
                    #print(j)
                    xerr.append(rsuf[j])
            else:
                xerr.append(serr[i])
        #print(sol[i],serr[i])
    x.append(1.0) #error scale
    xerr.append(0.1)
    if modeltype[3]==1:
        x.append(favg) #ascale
        xerr.append(0.1)
        x.append(1.0) #lscale
        xerr.append(0.1)
    x=np.array(x)
    xerr=np.array(xerr)

    ndim = len(x)

    p0=np.zeros((nwalkers,len(x)))
    k=-1
    for j in range(npar_sys):
        if np.abs(serr[j]) > 1.0e-30:
            k=k+1
            for ii in range(nwalkers):
                maxiter=10
                iter=0
                p0[ii,k]=1.0e30
                while p0[ii,k] >= rmax[j]+eps or p0[ii,k] <= rmin[j]-eps:
                    iter=iter+1
                    if iter < maxiter:
                        p0[ii,k]=x[k]+np.random.normal(scale=xerr[k])
                    else:
                        p0[ii,k]=np.min([np.max([x[k]+np.random.normal(scale=xerr[k]),rmin[j]]),\
                            rmax[j]])   
    for i in range(nbodies):
        for j in range(npar_sys):
            n=npar_sys+i*npar_pln+j
            if np.abs(serr[n]) > 1.0e-30:
                k=k+1
                for ii in range(nwalkers):
                    maxiter=10
                    iter=0
                    p0[ii,k]=1.0e30
                    while p0[ii,k] >= rmax[j+npar_sys]-eps or p0[ii,k] <= rmin[j+npar_sys]+eps:
                        iter=iter+1
                        if iter < maxiter:
                            p0[ii,k]=x[k]+np.random.normal(scale=xerr[k])
                        else:
                            p0[ii,k]=np.min([np.max([x[k]+\
                                                     np.random.normal(scale=xerr[k]),\
                                                     rmin[j+npar_sys]]),rmax[j+npar_sys]])
                            #print(ii,k,p0[ii,k],x[k],rmin[j+npar_sys],rmax[j+npar_sys])
    
    if modeltype[3]==0:
        p0[:,[len(x)-1]]=[np.random.rand(1)+0.5 for ii in range(nwalkers)] #deal with dscale
    elif modeltype[3]==1:
        p0[:,[len(x)-3]]=[np.random.rand(1)+0.5 for ii in range(nwalkers)] #deal with dscale
        p0[:,[len(x)-2]]=[(np.random.rand(1)+0.5)*favg for ii in range(nwalkers)] #deal with ascale
        obslen,mediandt=get_dtime(phot)
        tdur=tfit5.transitdur(sol,1)/86400.0 #Transit duration in days.
        p0[:,[len(x)-1]]=[np.random.rand(1)*(2*tdur-2*mediandt)+2*mediandt \
          for ii in range(nwalkers)] #deal with lscale
    p0[0,:]=np.copy(x) #make our good solution a walker
    
    return x,p0,ndim

def setup_priors(sol,nhyper):
    """Setup default priors for the transit model.
    Note: rmin and rmax contain parameter boundaries for valid models."""

    priors=[None]*(len(sol)+int(nhyper))

    npar_sys=8 #number of system-wide parameters (e.g., limb-darkening)
    npar_pln=10 #planet specific parameters (e.g., r/R*)

    npl=int(len(sol)) #number of parameters
    for i in range(npl):
        if (np.mod(i-npar_sys,npar_pln)==0) & (i>=npar_sys):
            priors[i]=[sol[i],10.0,10.0] #generous prior on EPO
        if (np.mod(i-npar_sys,npar_pln)==1) & (i>=npar_sys):
            priors[i]=[sol[i],10.0,10.0] #generous prior on PER


    return priors

def addeccn(serr):

    npar_sys=8 #number of system-wide parameters (e.g., limb-darkening)
    npar_pln=10 #planet specific parameters (e.g., r/R*)

    npl=int(len(serr)) #number of parameters
    for i in range(npl):
        if (np.mod(i-npar_sys,npar_pln)==4) & (i>=npar_sys):
            serr[i]=0.003
        if (np.mod(i-npar_sys,npar_pln)==5) & (i>=npar_sys):
            serr[i]=0.003

    return serr

def calcdiffs(phot1,phot2):

    npt1=len(phot1.time)
    npt2=len(phot2.time)
    diffs=np.zeros((npt1,npt2))
    for i in range(npt1):
        for j in range(npt2):
            diffs[i,j]=np.abs(phot1.time[i]-phot2.time[j])

    return diffs

def setup_gibbs():
    #setting up walkers.  First array gives amount of shuffle to add to good solution.  Second is min bounds on
    #variable, third is max bounds on variable.
    #               rho.     c1     c2.     q1    q2.    dil   vof   zpt.     epo.   per    b    rdr     ecw    esw
    #              krv     ted    ell    alb    DSC,     a       l
    rsuf=np.array([ 0.003 ,3.0e-4,3.0e-4,3.0e-4,3.0e-4,0.001, 0.1,   1.0e-6, 1.0e-6,1.0e-8, 0.01,1.0e-6, 0.001, 0.001,\
                    0.1,    0.1,  0.1,   0.1,   0.5,    1.0,    1.0])
    rmin=np.array([ 1.0e-5,-10.0 ,-10.0 ,0.0000,0.0000,0.000,-1.0e6,-1.0   ,-1.0e+5,1.0e-6, 0.0 ,0.0   ,-1.000,-1.000,\
                   -1.0e6,-1.0e6,-1.0e6,-1.0e6, 0.0   , 0.0,    0.0])
    rmax=np.array([ 1.0e3 , 10.0 , 10.0 ,1.0000,1.0000,1.000, 1.0e6, 1.0   , 1.0e+5,1.0e+6,10.0 ,2.0   , 1.000, 1.000,\
                    1.0e6, 1.0e6, 1.0e6, 1.0e6, 1.0e10, 1.0e10, 1.0e10 ])

    #DSC - scale for photometry errors
    #l - length scale for GPs

    return rsuf,rmin,rmax;


def transitmodel (sol,time, itime=-1, ntt=0, tobs=0, omc=0, dtype=0 ): 
    "read in transitmodel solution"  
    nplanet=int((len(sol)-8)/10) #number of planets
    
    if type(itime) is int :
        if itime < 0 :
            itime=ones(len(time))*0.020434
        else:
            itime=ones(len(time))*float(itime)
    
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

    tmodel= np.zeros(len(time)) #contains the transit model
    tfit5.transitmodel(nplanet,sol,time,itime,nttin,tobsin,omcin,tmodel,dtypein)
    return tmodel;


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

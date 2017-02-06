#!/usr/bin/env python

import sys,os, string
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting, polynomial

import pyfits
import scipy

import time
start_time1 = time.time()

#Splitting the looong file into smaller files
#
#from itertools import izip_longest
#
#def grouper(n, iterable, fillvalue=None):
#    "Collect data into fixed-length chunks or blocks"
#    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
#    args = [iter(iterable)] * n
#    return izip_longest(fillvalue=fillvalue, *args)
#
#n = 13323
#
#with open('test4.hires1b') as f:
#    for i, g in enumerate(grouper(n, f, fillvalue=''), 1):
#        with open('small_file_{0}'.format(i * n), 'w') as fout:
#            fout.writelines(g)
#

#reading the file


SED_FILE=[]
for SED_FILE_t in os.listdir('.'):
    if SED_FILE_t.startswith("small_file"):
        #print(SED_FILE_t)
        SED_FILE.append(SED_FILE_t)



####
#
#   Funcion lee_spectra
#   lee un spectrum y lo deja como x e y
#
#
#
#####

def Lee_spectra(spectrum):
    #lee un spectrum y lo deja como x e y

    imagen=pyfits.getdata(spectrum,header=False)
    header=pyfits.getheader(spectrum)

    #empty array
    Lambda_t=[]
    Flux_t=[]

    for i in range (len(imagen)):
        y=imagen[i]
        x=i*header['CDELT1']+header['CRVAL1']
        Lambda_t.append(float(imagen[i]))
        Flux_t.append(float(i*header['CDELT1']+header['CRVAL1']))
        #print x,y

    Flux=np.array(Lambda_t)
    Lambda=np.array(Flux_t)
    x2=[]
    y2=[]
    k=0
    for j in range(len(Lambda)):
        if (Lambda[j] > 3816) and (Lambda[j] <= 94260):
            x2.append(Lambda[j])
            y2.append(Flux[j])
            k=k+1
        else:
            j=j+1
        
    return x2,y2

    #x=Lambda
    #y=Flux
    #return x,y




#######
# poly_fit, fitea un polinomio a datos
# xp e yp correspondend a x e y a ser fiteado
#
# en este caso corresponden al x e y del output de la
# region discontinua
#
#
#######

def poly_fit(xp,yp,grado_pol):
    
    t_init = polynomial.Polynomial1D(degree=int(grado_pol))
    fit_t = fitting.LevMarLSQFitter()
    t = fit_t(t_init, xp, yp)
    return t

#########################



def lee_archivo(archivo,minimo,maximo):
    TIME=[]
    WAVE=[]
    LOG_LUM=[]
    FNORM=[]
    for line in open(archivo):
        li=line.strip()
        if not li.startswith("#"):##

            todo=line.split()
            TIME.append(float(todo[0]))
            WAVE.append(float(todo[1]))
            LOG_LUM.append(float(todo[2]))
            FNORM.append(float(todo[3]))
            # data is sorted in order to do the KDE. and [np.newaxis] is added.
            # The KDE needs this format (sight!)
            #val2=np.array(np.sort(value))[np.newaxis]
            #val=val2.transpose()
            yr_t=np.array(TIME) 
            wav_t=np.array(WAVE)+4.558#+4.558 
            lum_t=np.array(LOG_LUM)
            flux_t=np.array(FNORM)
    yr_ta=[]
    wav_ta=[]
    lum_ta=[]
    flux_ta=[]
    #print len(yr_t)
    k=0
    for i in range (len(yr_t)):
        if (wav_t[i] >= minimo-1) and (wav_t[i]<=maximo+1):         
            
            yr_ta.append(yr_t[i]) 
            wav_ta.append(wav_t[i]) 
            lum_ta.append(lum_t[i]) 
            flux_ta.append(flux_t[i])

            yr_tb  =np.array(yr_ta  )
            wav_tb =np.array(wav_ta )
            lum_tb =np.array(lum_ta )
            flux_tb=np.array(flux_ta)

            k=k+1
        else:
            i=i+1
    return(yr_tb,wav_tb,lum_tb,flux_tb)



x_sci,y_sci=Lee_spectra(sys.argv[1])

#Minimization a la Brute force!


Minimo=9E999
for l in range(len(SED_FILE)):

    yr,wav,lum,flux=lee_archivo(SED_FILE[l],min(x_sci),max(x_sci))
    #fitting a line
    Y_FIT=poly_fit(wav,flux,1)

    #interpolating



    flux_N=(flux/Y_FIT(wav))
    inter_low=np.interp(x_sci,wav,flux_N)
    inter_hi=np.interp(wav,x_sci,y_sci)


    #inter_low=np.interp(x_sci,wav,flux)
    #inter_hi=np.interp(wav,x_sci,y_sci)
    #y_resta_hires=inter_hi-wav



    y_resta_lores=-inter_low+y_sci
    y_minimizar= sum(np.sqrt(y_resta_lores*y_resta_lores))
    chi2=sum((-inter_low+y_sci)*(-inter_low+y_sci)/inter_low)

    if (y_minimizar < Minimo):
        Minimo=y_minimizar
        Name=SED_FILE[l]
        print Minimo,SED_FILE[l],chi2,yr[0]
    else:
        Minimo=Minimo
        print Minimo,SED_FILE[l],y_minimizar,chi2,yr[0]

    file_plot=str(sys.argv[1][:-10])+'_'+str(yr[0])[:-8]+'.pdf'
    plt.plot(x_sci,y_sci)
    plt.plot(x_sci,inter_low,color='r')
    plt.plot(x_sci,y_resta_lores)
    #file_plot=str(sys.argv[1][:-10])+'_'+str(yr[0][:-8])+'.pdf'
    
    plt.xlabel('$\lambda$ [Angstroms]')
    plt.ylabel('Normalized')
    plt.savefig(file_plot)
    plt.close()
    #file_plot
#print len(inter_hi)


#print data2[0],test


#plt.plt(wav,RESTA)

#plt.plot(data1[0],test)
#plt.plot(x_sci,test)
#plt.plot(x_sci,y_sci)
#plt.plot(x_sci,y_resta_lores)
print "minimum value SED file name is\n"
print Name
print "Total execution time: %.2f min" % ((time.time() - start_time1)/60)

#plt.show()

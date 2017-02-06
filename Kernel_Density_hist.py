#!/usr/bin/env python

#Read a single column file, creates a normalized histogrtam and fit a probability density function (Kernel Densite Estimator)

import sys,os,string,time
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting, polynomial




def lee_archivo(archivo):
    value=[]
    for line in open(archivo):
        li=line.strip()
        if not li.startswith("#"):

            todo=line.split()
            value.append(float(todo[0]))
            # data is sorted in order to do the KDE. and [np.newaxis] is added.
            # The KDE needs this format (sight!)
            val2=np.array(np.sort(value))[np.newaxis]
            val=val2.transpose()
    return(val)


vel=lee_archivo("SL538_900sORDER89_END_voigt_bootN.txt") #"File Name"

#vel=vel2.reshape((1,-1))
from sklearn.neighbors.kde import KernelDensity
#kde = KernelDensity(kernel='gaussian', bandwidth=0.0001,leaf_size=60).fit(vel)
kde2 = KernelDensity(kernel='gaussian', bandwidth=0.005).fit(vel)
y_kde=kde.score_samples(vel)
y_kde2=(kde2.score_samples(vel))
#y_kde2=(kde2.score(vel))



#print vel
print kde
print kde.score_samples(vel)


a,b,c=plt.hist((vel),90,normed=1)#xs=np.linspace(0,8,200)
 
#print a,b,c
plt.plot(vel,np.exp(y_kde2))#*(np.log(len(vel))))

plt.draw()
#time.sleep(0.05)
plt.show()

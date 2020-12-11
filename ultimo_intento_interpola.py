import pandas as pd
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from operator import attrgetter
from sklearn.preprocessing import normalize
from scipy.interpolate import interp1d,interpn,interp2d,Rbf
from scipy.interpolate import griddata



def todo():
    global interpola
    #No interpolation
    #age= [1.0,2.0,3.0,4.0,5.0,6.0] #NO interpolation on AGE #se ordena por este factor
    geom_fact=[0.03,3]             #ordednar 
    cut=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]   #Ordenar

    #para interpolar
    
    log_U=[-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0]                  #extraporlar en 150+1
    OH12= [6.6,6.8,7.0,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.0,9.2,9.4] #extrapolar 28+1. la idea es que quede con 0.1
    NOlog=[-2, -1.75,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25 , 0.0]    #extrapolar por 40+1 dejarlo cada 0.05
    #print(len(NOlog))
    #print(NOlog[1])
    #Mostrar todo

    data=pd.read_csv('BOND.csv',header=0, names= ['AGE','geom_factor','lU_mean','nebula_cut','OH','NO','logU','O2','O3_4363','O3','N2','S2'] )
    df=pd.DataFrame(data,columns= ['AGE','geom_factor','lU_mean','nebula_cut','OH','NO','logU','O2','O3_4363','O3','N2','S2'] )#
    df=df.round({'nebula_cut':1})

    dt=pd.DataFrame(data,columns= ['AGE','geom_factor','lU_mean','nebula_cut','OH','NO','logU','O2','O3_4363','O3','N2','S2'] )#
    dt=dt.round({'nebula_cut':1})


    age=[1000000.0,2000000.0,3000000.0,4000000.0,5000000.0,6000000.0] #NO interpolation on AGE #se ordena por este factor
    geom_fact=[0.03,3]             #ordednar 
    cut=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]   #Ordenar
    log_U=[-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0]                  #extraporlar en 150+1
    

    i_lum=np.round(np.linspace(np.min((df['lU_mean'])),np.max((df['lU_mean'])),num=101,endpoint=True),decimals=2)
    i_oh=(np.linspace(np.min((df['OH'])),np.max((df['OH'])),num=57,endpoint=True))
    i_no=(np.linspace(np.min((df['NO'])),np.max((df['NO'])),num=41,endpoint=True))
    delta_no=i_no[1]-i_no[0]
    delta_oh=i_oh[1]-i_oh[0]
    delta_lum=i_lum[1]-i_lum[0]
    
    #grid 
    x,y=np.mgrid[np.min(i_oh):np.max(i_oh):delta_oh,np.min(i_no):np.max(i_no):delta_no]
    
    
    for a in range(len(age)):
        for g in range(len(geom_fact)):
            for c in range(len(cut)):
                for l in range(len(log_U)): 
                     A=df.loc[(df['geom_factor']==geom_fact[g])  & (df['AGE']==age[a]) & (df['nebula_cut']==cut[c]) &(df['lU_mean']==log_U[l])]# 

                     #Estos son los puntos a extrapolar
                     #x,y,z=np.mgrid[np.min(i_oh):np.max(i_oh):delta_oh,np.min(i_no):np.max(i_no):delta_no,np.min(i_lum):np.max(i_lum):delta_lum]

                     f1=np.round(griddata((A['OH'],A['NO']),A['AGE'], (x,y),method='cubic'),decimals=1)   
                     f2=np.round(griddata((A['OH'],A['NO']),A['geom_factor'], (x,y),method='cubic'),decimals=2)   
                     f3=np.round(griddata((A['OH'],A['NO']),A['lU_mean'], (x,y),method='cubic'),decimals=2)   
                     f4=np.round(griddata((A['OH'],A['NO']),A['nebula_cut'], (x,y),method='cubic'),decimals=1)   
                     f5=np.round(griddata((A['OH'],A['NO']),A['OH'], (x,y),method='cubic'),decimals=2)   
                     f6=np.round(griddata((A['OH'],A['NO']),A['NO'], (x,y),method='cubic'),decimals=2)   
                     f7=np.round(griddata((A['OH'],A['NO']),A['logU'], (x,y),method='cubic'),decimals=20)   
                     f8=np.round(griddata((A['OH'],A['NO']),A['O2'], (x,y),method='cubic'),decimals=20)   
                     f9=np.round(griddata((A['OH'],A['NO']),A['O3_4363'], (x,y),method='cubic'),decimals=20)   
                     f10=np.round(griddata((A['OH'],A['NO']),A['O3'], (x,y),method='cubic'),decimals=20)   
                     f11=np.round(griddata((A['OH'],A['NO']),A['N2'], (x,y),method='cubic'),decimals=20)   
                     f12=np.round(griddata((A['OH'],A['NO']),A['S2'], (x,y),method='cubic'),decimals=20)   
                     #print(len(x),len(y),f10.shape)
                     l1,l2=f10.shape
                     for i in range(l1):
                         for j in range(l2):
                             print (f1[i,j],f2[i,j],f3[i,j],f4[i,j],f5[i,j],f6[i,j],f7[i,j],f8[i,j],f9[i,j],f10[i,j],f11[i,j],f12[i,j])

                   

    return()

todo()

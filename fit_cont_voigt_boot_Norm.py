#!/User/bin/env python

### Idea: 
#- hacer bootstraping con los fiteos de la gaussiana
#- Tener los errores
#- plot con los intervalos

#importing system commands
import sys,os,string

#scientific packages
import pyfits
import scipy
from scipy.optimize import curve_fit, leastsq
import numpy as np
from numpy import random, exp,sqrt
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting, polynomial

N_ITER=15005 # Numero de iteraciones

imagen_in=sys.argv[1]
l_min_izq=float(sys.argv[2])
l_max_izq=float(sys.argv[3])

l_min_der=float(sys.argv[4])
l_max_der=float(sys.argv[5])

guess_line=float(sys.argv[6])
guess_FWHM=float(sys.argv[7])


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

    imagen=pyfits.getdata(imagen_in,header=False)
    header=pyfits.getheader(imagen_in)

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
    x=Lambda
    y=Flux
    return x,y


x_sci,y_sci=Lee_spectra(imagen_in)

#########################
#
# Funcion Region
# Toma una region de un espectro entre
# un minimo lambda y un maximo lambda
# x e y corresponden a lamba y cuentas o flujo
#
##########################

def region(minimo,maximo,x,y):
        xar=[]
        yar=[]
        for i in range(len(x)):
            if (x[i] > minimo) and (x[i] <maximo):
                xar.append(float(x[i]))
                yar.append(float(y[i]))

        xar=np.array(xar)
        yar=np.array(yar)
        return xar,yar


#########################
#
# Funcion Region_discontinuo
# Toma dos regiones de un espectro separadoas entre
# un minimo lambda y un maximo lambda a la izquierda y un 
# un minimo lambda y un maximo lambda a la deracha de la emission o obsorption 
# x e y corresponden a lamba y cuentas (o flujo)
#
##########################

def region_discontinua(minimo1,maximo1,minimo2,maximo2,x,y):
    xar=[]
    yar=[]
    for i in range(len(x)):
        if ((x[i] > minimo1) and (x[i] <maximo1)) or ((x[i] > minimo2) and (x[i] <maximo2)):
            xar.append(float(x[i]))
            yar.append(float(y[i]))
            
    xar=np.array(xar)
    yar=np.array(yar)
    return xar,yar

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

#calclulo original

x=x_sci 
y=y_sci 

#xa1,ya1=region(l_min_izq,l_max_izq,x,y)
#xa2,ya2=region(l_min_der,l_max_der,x,y)

############

xspec_o,yspec_o=region(l_min_izq,l_max_der,x,y)

#################


###Fitting regions with a polynomio


x_cont_o,y_cont_o=region_discontinua(l_min_izq,l_max_izq,l_min_der,l_max_der,x,y)


#cont1=poly_fit(xa1,ya1,12)
#cont2=poly_fit(xa2,ya2,12)
cont3_o=poly_fit(x_cont_o,y_cont_o,12)

#print cont1
#print cont2

#res1= -cont1(xa1)+ ya1
#res2= -cont2(xa2)+ ya2
res3_o= -cont3_o(x_cont_o)+ y_cont_o

#se aplica el polinomio al espectro en la zona de interes
res4_o= -cont3_o(xspec_o)+yspec_o

#####################
#
# Normalization!!!
#
######################
res4_o= yspec_o/cont3_o(xspec_o)


######################


################

#iteracion 1
t_init4_o = models.GaussianAbsorption1D(amplitude=1, mean=guess_line, stddev=guess_FWHM)
fit_t4_o = fitting.LevMarLSQFitter()
t4_o = fit_t4_o(t_init4_o, xspec_o, res4_o)

a_science=t4_o.mean.value
b_science=t4_o.stddev.value
Amplitud=t4_o.amplitude.value

#iteracion 2 para gaussiana
t_init4_o = models.GaussianAbsorption1D(amplitude=Amplitud, mean=a_science, stddev=b_science)
fit_t4_o = fitting.LevMarLSQFitter()
t4_o = fit_t4_o(t_init4_o, xspec_o, res4_o)

residuo_o=-t4_o(xspec_o)+res4_o

#print "resultados",t_init4,t4
print "resultados",t4_o
a_science=t4_o.mean.value
b_science=t4_o.stddev.value
c_science_amplitude=t4_o.amplitude.value

#Aplicando FWHM guess
guess_FWHM_gauss=-float(-b_science)
print guess_FWHM_gauss
#exit(0)


##print t4.mean
Lambda_gauss_fit_sci="{:10.3f}".format(a_science)
Sigma_gauss_fit_sci="{:10.3f}".format(b_science)
print  Lambda_gauss_fit_sci, Sigma_gauss_fit_sci
#print a ,b



LINE=[]
FWHM=[]
aa=[]
bb=[]
guess_line=t4_o.mean.value

############################################
#
#
#   Definiciones de Voigt 
#
#
############################################
from matplotlib.pyplot import figure, show, rc
from scipy.special import wofz
from kapteyn import kmpfit

ln2 = np.log(2)

def voigt(x, y):
   # The Voigt function is also the real part of 
   # w(z) = exp(-z^2) erfc(iz), the complex probability function,
   # which is also known as the Faddeeva function. Scipy has 
   # implemented this function under the name wofz()
   z = x + 1j*y
   I = wofz(z).real
   return I

def Voigt(nu, alphaD, alphaL, nu_0, A, a=0, b=0):
   # The Voigt line shape in terms of its physical parameters
   f = np.sqrt(ln2)
   x = (nu-nu_0)/alphaD * f
   y = alphaL/alphaD * f
   backg = a + b*nu 
   V = A*f/(alphaD*np.sqrt(np.pi)) * voigt(x, y) + backg
   return V

def funcV(p, x):
    # Compose the Voigt line-shape
    alphaD, alphaL, nu_0, I, a, b = p
    return Voigt(x, alphaD, alphaL, nu_0, I, a, b)

def funcV2(p, x):
    # Compose the Voigt line-shape
    alphaD, alphaL, nu_0, I, a, b = p
    return Voigt(x, alphaD, alphaL, nu_0, I, a, b)

def residualsV(p, data):
   # Return weighted residuals of Voigt
   x, y, err = data
   return (y-funcV(p,x)) / err


N = len(res4_o)
err = np.ones(N)
A = t4_o.amplitude.value
alphaD = 0.5
alphaL = 0.5
a = 6
b = 0
nu_0 = guess_line
p0 = [alphaD, alphaL, nu_0, A, a, b]



fitter = kmpfit.Fitter(residualsV, data=(xspec_o,res4_o,err))
fitter.parinfo = [{}, {}, {}, {}, {}, {'fixed':True}]  # Take zero level fixed in fit
fitter.fit(params0=p0)
#print "\n========= Fit results Voigt profile =========="
#print "Initial params:", fitter.params0
#print "Params:        ", fitter.params
#print "Iterations:    ", fitter.niter
#print "Function ev:   ", fitter.nfev 
#print "Uncertainties: ", fitter.xerror
#print "dof:           ", fitter.dof
#print "chi^2, rchi2:  ", fitter.chi2_min, fitter.rchi2_min
#print "stderr:        ", fitter.stderr   
#print "Status:        ", fitter.status
##empirical values for voight:
c1 = 1.0692
c2 = 0.86639
hwhm = 0.5*(c1*alphaL+np.sqrt(c2*alphaL**2+4*alphaD**2))
print "\nFWHM Voigt profile:     ", 2*hwhm

#print fitter
#print fitter.message
print "\n Central lambda of the feature: ",fitter.params[2]


############################
#
#segunda iteracion Voight!
#
###############################

N = len(res4_o)
err = np.ones(N)
A = float(fitter.params[3])#-2000
alphaD = float(fitter.params[0])#0.5
alphaL = float(fitter.params[1])#0.5
a = float(fitter.params[4])#6
b = 0
nu_0 = fitter.params[2]
p0 = [alphaD, alphaL, nu_0, A, a, b]



fitter = kmpfit.Fitter(residualsV, data=(xspec_o,res4_o,err))
fitter.parinfo = [{}, {}, {}, {}, {}, {'fixed':True}]  # Take zero level fixed in fit
fitter.fit(params0=p0)
#print "\n========= Fit results Voigt profile =========="
#print "Initial params:", fitter.params0
#print "Params:        ", fitter.params
#print "Iterations:    ", fitter.niter
#print "Function ev:   ", fitter.nfev 
#print "Uncertainties: ", fitter.xerror
#print "dof:           ", fitter.dof
#print "chi^2, rchi2:  ", fitter.chi2_min, fitter.rchi2_min
#print "stderr:        ", fitter.stderr   
#print "Status:        ", fitter.status
#empirical values for voight:
c1 = 1.0692
c2 = 0.86639
hwhm = 0.5*(c1*alphaL+np.sqrt(c2*alphaL**2+4*alphaD**2))
print "\nFWHM Voigt profile:     ", 2*hwhm

#print fitter
#print fitter.message
print "\n Central lambda of the feature iter2: ",fitter.params[2]

#parametros=fitter.params

############################
#
#TERCERA iteracion Voight!
#
###############################

N = len(res4_o)
err = np.ones(N)
A = float(fitter.params[3])#-2000
alphaD = float(fitter.params[0])#0.5
alphaL = float(fitter.params[1])#0.5
a = float(fitter.params[4])#6
b = 0
nu_0 = fitter.params[2]
p0 = [alphaD, alphaL, nu_0, A, a, b]



fitter = kmpfit.Fitter(residualsV, data=(xspec_o,res4_o,err))
fitter.parinfo = [{}, {}, {}, {}, {}, {'fixed':True}]  # Take zero level fixed in fit
fitter.fit(params0=p0)
#print "\n========= Fit results Voigt profile =========="
#print "Initial params:", fitter.params0
#print "Params:        ", fitter.params
#print "Iterations:    ", fitter.niter
#print "Function ev:   ", fitter.nfev 
#print "Uncertainties: ", fitter.xerror
#print "dof:           ", fitter.dof
#print "chi^2, rchi2:  ", fitter.chi2_min, fitter.rchi2_min
#print "stderr:        ", fitter.stderr   
#print "Status:        ", fitter.status
#empirical values for voight:
c1 = 1.0692
c2 = 0.86639
hwhm = 0.5*(c1*alphaL+np.sqrt(c2*alphaL**2+4*alphaD**2))
print "\nFWHM Voigt profile:     ", 2*hwhm

#print fitter
#print fitter.message
print "\n Central lambda of the feature iter3: ",fitter.params[2]

parametros=fitter.params
#antes del random 



#fin 




##Aplicando Random 


for N in range(N_ITER):
    x=x_sci 
    y=y_sci +  np.sqrt((y_sci))*random.randn(len(x_sci))
    #print y[600], y_sci[600]
    #print "iter",N
    
    #xa1,ya1=region(l_min_izq,l_max_izq,x,y)
    #xa2,ya2=region(l_min_der,l_max_der,x,y)
    xspec,yspec=region(l_min_izq,l_max_der,x,y)

    ###Fitting regions with a polynomio


    x_cont,y_cont=region_discontinua(l_min_izq,l_max_izq,l_min_der,l_max_der,x,y)


    #cont1=poly_fit(xa1,ya1,12)
    #cont2=poly_fit(xa2,ya2,12)
    cont3=poly_fit(x_cont,y_cont,12)

    #print cont1
    #print cont2

    #res1= -cont1(xa1)+ ya1
    #res2= -cont2(xa2)+ ya2
    #################
    #
    #
    #
    #   Caombiado el 6 de Junio
    #
    #
    ########################

    #res3= -cont3(x_cont)+ y_cont
    #se aplica el polinomio al espectro en la zona de interes
    #res4= -cont3(xspec) + yspec


    res3= y_cont/cont3(x_cont)

    #se aplica el polinomio al espectro en la zona de interes
    res4= yspec/cont3(xspec)




    ################
        

    fitterN = kmpfit.Fitter(residualsV, data=(xspec,res4,err))
    fitterN.parinfo = [{}, {}, {}, {}, {}, {'fixed':True}]  # Take zero level fixed in fit
    fitterN.fit(params0=p0)
    #print "\n========= Fit results Voigt profile =========="
    #print "Initial params:", fitter.params0
    #print "Params:        ", fitter.params
    #print "Iterations:    ", fitter.niter
    #print "Function ev:   ", fitter.nfev 
    #print "Uncertainties: ", fitter.xerror
    #print "dof:           ", fitter.dof
    #print "chi^2, rchi2:  ", fitter.chi2_min, fitter.rchi2_min
    #print "stderr:        ", fitter.stderr   
    #print "Status:        ", fitter.status


    #print fitter.params
    #print fitter.message


    ##seleccionando la region que se fitea solamente
    
    #residuo=funcV(fitterN.params,xspec)-res4
    residuo=funcV(fitterN.params,xspec)-res4
    residuo_x,residuo_y=region(l_max_izq,l_min_der,xspec,residuo)



    #t_init4 = models.GaussianAbsorption1D(amplitude=1, mean=guess_line, stddev=guess_FWHM)
    #fit_t4 = fitting.LevMarLSQFitter()
    #t4 = fit_t4(t_init4, xspec, res4)

    #residuo=-t4(xspec)+res4

    #print "resultados",t_init4,t4
    #print "here",fitter.params[2]
    aa.append(float(fitterN.params[2]))
    #b.append(float(t4.stddev.value))
    
    ##print t4.mean
    Lambda_fit="{:10.3f}".format(aa[N])
    #Sigma_fit="{:10.3f}".format(b[N])
    #print  Lambda_fit, Sigma_fit
    #print a ,b
    
    LINE=np.array(aa)
   # FWHM=np.array(b)

name=imagen_in[:-21]

minimo_lambda=min(LINE)
maximo_lambda=max(LINE)
#min_sigma=min(FWHM)
#max_sigma=float(max(FWHM))

#a_science
#b_science

#ACA!!!
print "resultados",t4_o

#print fitter.params[2]
a_Moffat_science=fitter.params[2]

print fitter.params
#b_science=t4_o.stddev.value
#c_science_amplitude=t4_o.amplitude.value

##print t4.mean
Lambda_Moffat_fit_sci="{:10.3f}".format(a_Moffat_science)
#Sigma_fit_sci="{:10.3f}".format(b_science)
print  Lambda_Moffat_fit_sci
#print a ,b





e_min=( maximo_lambda-a_science)
e_max=(-minimo_lambda+a_science)
#e_min=float(max_sigma)-float(Lambda_fit_sci)

File_out2=imagen_in[:-5]+'_voigt_boot_TEST_histogram.png'
File_text=imagen_in[:-5]+'_voigt_boot.txt'
np.savetxt(File_text,LINE,fmt='%.5f')
File_result=imagen_in[:-5]+'_voigt_boot_result.txt'

f1=open(File_result,'w')
#f1.write(Lambda_fit_sci)
f1.write(str(Lambda_Moffat_fit_sci)+' '+str(minimo_lambda)+' '+str(maximo_lambda)+'\n ')
#f1.write(maximo_lambda)
f1.write(str(Lambda_Moffat_fit_sci) +' +'+str(e_min) +' -'+str(e_max)+'\n') 
f1.close() 

print  Lambda_Moffat_fit_sci#, Sigma_fit_sci
print  Lambda_Moffat_fit_sci,"min=",minimo_lambda,"max=",maximo_lambda 
#print  Sigma_fit_sci,"min=",min_sigma,"max=",max_sigma 
print  Lambda_Moffat_fit_sci,"+",e_min,"-",e_max 




#bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
#plt.hist(LINE,10)
plt.hist(LINE,bins=10)
#plt.plot(n,bines)
plt.xlabel('Angstroms')
plt.ylabel('N')
plt.savefig(File_out2)
plt.close() 




plt.subplot(311)
plt.title(name) 
plt.ylabel('Conuts')
plt.xlim(l_min_izq,l_max_der)
plt.plot(xspec_o,yspec_o,'k-',lw=1, label='')
#plt.plot(xa1, ya1, 'b-', lw=2, label='Trapezoid')
#plt.plot(xa1,cont1(xa1))
#plt.plot(xa2, ya2, 'r-', lw=2, label='Gaussian')
#plt.plot(xa2,cont2(xa2))
#plt.plot(x_cont,cont3(x_cont))
plt.plot(xspec_o,cont3_o(xspec_o),'c-',lw=1, label='')
###


plt.subplot(312)
#plt.plot(xa1, res1, 'b-', lw=2, label='Trapezoid')
#plt.plot(xa2, res2, 'r-', lw=2, label='Gaussian')
#plt.plot(x_cont, res3, 'g-', lw=2, label='Gaussian')
plt.xlim(l_min_izq,l_max_der)
plt.ylabel('Difference')
plt.plot(xspec_o,res4_o, 'c-',lw=1,label='gaus')
plt.plot(xspec_o,funcV(fitter.params,xspec),'g-', lw=2,label='gauss')

parametros=fitter.params

#residuoV=funcV(fitter.params,xspec_o)-res4
#residuoV=funcV(parametros,xspec_o)-res4_o
residuoV=res4_o/funcV(parametros,xspec_o)

#residuo_x,residuo_y=region(l_max_izq,l_min_der,xspec_o,residuoV)
residuo_x,residuo_y=region(l_min_izq,l_max_der,xspec_o,residuoV)



e_min_str="{:10.3f}".format(-e_min)
e_max_str="{:10.3f}".format(+e_max)

#plt.annotate(r'$\lambda$='+Lambda_Moffat_fit_sci+''+str(e_min_str)+''+str(e_max_str),xy=(l_min_izq,-1000))
plt.annotate(r'$\lambda$='+Lambda_Moffat_fit_sci+''+str(e_min_str)+''+str(e_max_str),xy=(l_min_izq,0.9))

#plt.plot(xspec,residuo,'g-', lw=2,label='gauss')
#plt.plot(xa2,cont2(xa2))
#


plt.subplot(313)
plt.xlim(l_min_izq,l_max_der)
plt.plot(xspec_o,residuoV,'g-', lw=2,label='gauss')
plt.ylabel('Residual')
plt.xlabel('$\lambda$ [Angstroms]')
#plt.plot(xa2,cont2(xa2))

#removing the .fits
File_out=imagen_in[:-5]+'_gaussfit_boot_voigt.png'

plt.savefig(File_out)

#plt.savefig('temp_'+str(int(N))+'.png')
#plt.show()
plt.close() #asi se borran los graficos!

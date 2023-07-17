# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:59:13 2020

@author: VeroGeorl
"""
import numpy as np
import scipy.signal as sc
import matplotlib.pyplot as plt
#import plotly.plotly as py
#import plotly.graph_objs as go
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from math import* 
from scipy import *
import fAttenuations_Gl as f_att

""" Paramètres de simulation """ 
C=3e8*1e-9 #vitesse de la lumière et en nanoseconde (m/ns)
rho= 0.8 #coefficient de reflexion
Ts=20#gain du filtre 
n=1.5#indice de réfraction



"""-----Emetteur----"""
 
theta = 60*(np.pi/180)  #demi angle en degré 
m=-np.log10(2)/np.log10(np.cos(theta)) #Mode lambertien 
P_Total=2000 #puissance totale transmise en mw


"""----Recepteur----"""
Adet=np.exp(-4)#aire de detection du recepteur
FOV=60*(np.pi/180) #field of view
g=np.square(n)/np.square(np.sin(FOV))#gain 
R=0.65

"""----Dimensionnement de la piece ----"""
lx=3
ly=1.5
lz=2.1
KK=10
#dimension de la pièce
Nx=int(lx*KK)
Ny=int(np.round(ly*KK))
Nz=np.round(lz*KK)
#nombre de grille pour chaque surface
dA=(lz*ly)/(Ny*Nz)
#aire de la grille

x_=np.linspace(-lx/2, lx/2, Nx)
y_=np.linspace(-ly/2, ly/2, Ny)
z=np.linspace(-lz/2, lz/2, int(Nz))

xr, yr=np.meshgrid(x_,y_)
#vi=0.01 #km
vi=0.01
lam=0.450 #nm

""" Attenuation """
A_fog,A_fog_dB=f_att.compute_attenuation_fog(vi,lam)
A_smoke,A_smoke_dB=f_att.compute_attenuation_smoke_v(vi)

"""Other parameters to compute SNR"""
e=1.6E-19 #C (A.s) Electron charge
kb=1.38E-23 #J/K Boltzmann constant
T=298 #K temperature

G=10 #openloop voltage gain
gm= 30 #mS transconductance
GG=1.5 #FET vhannel noise factor
BW_optF= 10 #nm BW of bandpass optical filter

bg_irra= 5.8 #µW/cm².nm


Capa=112 #picoF/cm² capacitance per unit area

I2= 0.562
I3= 0.868 # Noise bandwidth factors
D_r= 10E6 #Transmission data rate bit/s

"""----Simulation--------"""
TP1=np.array([0,0,lz/2])

H=(Nx, Ny)
h0=np.zeros(H)

delta_t=1/2
time_size=25
#t_vector=np.arange(time_size/delta_t)
#h_vector= np.zeros(len(t_vector))

H_VEC=np.zeros((Nx, Ny,int(time_size/delta_t))) #H_VEC is the matrix with all the impulse response


for ii in range(int(Nx)):
    
    for jj in range (int(Ny)):
        t_vector=np.arange(time_size/delta_t)
        h_vector= np.zeros(len(t_vector))
        RP = np.array([x_[ii], y_[jj], -lz/2])#position de Rx
        """-------- TRANSMISSION LOS ----------""" 
        D1=np.sqrt(np.vdot(TP1-RP, TP1-RP))#Vecteur distance entre Tx et Rx
        tau1=D1/C
        index1=np.argwhere(round(tau1/delta_t)==t_vector)
        cosphi=lz/D1
        #print(index1)
        if np.abs(np.arccos(cosphi))<=FOV:
            TEST00=h0[ii, jj]
            TEST01=((m+1)*Adet*(cosphi**(m+1)))/(2*np.pi*np.square(D1))
            h0[ii, jj]= TEST00 + TEST01
            P_rec_A1=h0*P_Total*Ts*g
            h_vector[index1]=h_vector[index1]+((m+1)*Adet*(cosphi**(m+1)))/(2*np.pi*np.square(D1))
        H_VEC[ii,jj]=np.asarray(h_vector)
            


"""This is the classical Plot of the power received on the floor"""
H=h0
P_rec=H*P_Total*Ts*g
P_rec_att=P_rec*np.exp(-(A_smoke+A_fog)*0.0025)

x_=np.linspace(0, lx, Nx)
y_=np.linspace(0, ly, Ny)
xr, yr=np.meshgrid(x_,y_)

fig = plt.figure(figsize=plt.figaspect(2.))
ax = Axes3D(fig)
im=ax.plot_surface(xr, yr, 10*np.log10(P_rec).transpose(), rstride=1, cstride=1, cmap=cm.plasma,linewidth=0, antialiased=False)
fig.colorbar(im,orientation="vertical")
#plt.axis('scaled')
plt.title('Level curves dBm') 
plt.figure()
min_=np.min(10*np.log10(P_rec))
max_=np.max(10*np.log10(P_rec))
vec=np.linspace(min_,max_,10)
cont=plt.contour(xr,yr,10*np.log10(P_rec).transpose(),vec,cmap=cm.plasma)
plt.clabel(cont)
plt.axis('scaled')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Power level curves at the bus shelter in dBm') 
fig.legend()
#print(H_VEC)

#%% 
"""This is the computation of the impulse response of the LOS contribution"""

#fig = plt.figure(figsize=plt.figaspect(2.))
fig = plt.figure(figsize=(7,5))

x=np.arange(0,time_size,delta_t)

#H=H_VEC[int((Nx-1)/2), int((Ny-1)/2)] #Where we check the impulse response
H=H_VEC[0,0] #This is the worst place to be
#H=H_VEC[15,7] #This is the worst place to be
RP = np.array([H[0], H[1], -lz/2])#position de Rx 

D=np.sqrt(np.vdot(TP1-RP, TP1-RP))*10**-3
num=H*np.exp(-(A_fog+A_smoke)*D)
num1=H*np.exp(-(A_smoke)*D)
num2=H*np.exp(-(A_fog)*D)
plt.stem(x,num*100/sum(H),markerfmt = 'bx', linefmt = 'g--', basefmt = 'm:',label='Smoke and Fog')
plt.stem(x,num1*100/sum(H),markerfmt = 'gv', linefmt = 'g--', basefmt = 'm:',label='Smoke')
plt.stem(x,num2*100/sum(H),linefmt='grey', markerfmt='D',basefmt = 'm:',label='Fog ')
plt.stem(x,H*100/sum(H),markerfmt = 'ro', linefmt = 'g--', basefmt = 'm:',label='Clear weather ')
max_clear=np.max(H*100/sum(H))
ind_clear=np.where(num == np.amax(num))
ind_clear=ind_clear[0]/2
max_smoke= np.max(num1*100/sum(H))
max_fog=np.max(num2*100/sum(H))
max_both=np.max(num*100/sum(H))
plt.annotate(str(np.round(max_smoke,2))+' %', xy=(ind_clear, max_smoke), xytext=(0, max_smoke),arrowprops=dict(facecolor='black', shrink=0.03), fontsize=16) #smoke
plt.annotate(str(np.round(max_both,2))+' %', xy=(ind_clear, max_both), xytext=(0, max_both),arrowprops=dict(facecolor='black', shrink=0.03), fontsize=16)
plt.annotate(str(max_clear)+' %', xy=(ind_clear, max_clear), xytext=(0,max_clear),arrowprops=dict(facecolor='black', shrink=0.03), fontsize=16)#clear weather
plt.annotate(str(np.round(max_fog,2))+' %', xy=(ind_clear, max_fog), xytext=(0, max_fog),arrowprops=dict(facecolor='black', shrink=0.03), fontsize=16) #fog
plt.xlabel('time [ns]',fontweight='bold', fontsize=16)
plt.title('Relative CIR for different weather (Visibility of '+ str(vi*1000)+ ' m)', fontsize=16) 
plt.xticks(size = 16)
plt.yticks(size = 16)
#plt.axes(fontsize=16)
plt.legend(fontsize=16)


#%%
""" Computation of the SNR curve """




"""
fig = plt.figure(figsize=plt.figaspect(2.))
RX1=np.array([lx/4,ly/4])
x=np.arange(0,25,0.5)
plt.stem(x,H_VEC[7,7]/sum(H_VEC[7,7]),'-',label='Impulse response center')
plt.stem(x,H_VEC[7,7]/sum(H_VEC[7,7]),'-',label='Impulse response center')
plt.xlabel('time in ns')
plt.legend()
plt.show()
"""
"""

x=np.arange(0,25,0.5)
plt.figure()
for i in range(len(H_VEC)):
    RP = np.array([x_[i], 0, -lz/2])#position de Rx
    D=np.sqrt(np.vdot(TP1-RP, TP1-RP))*10**-3
    num=H_VEC[i,0]*np.exp(-(A_fog+A_smoke)*D)
    #print(np.exp(-(A_fog+A_smoke)*D))
    plt.stem(x,num/sum(H_VEC[i,0]),'-',label='8, '+str(i))
    #plt.plot(x,num,'-',label='8, '+str(i))
    plt.xlabel('time in ns')
plt.legend()
plt.show()
"""
"""
fig = plt.figure(figsize=plt.figaspect(2.))
x=np.arange(0,25,0.5)
for i in range(len(H_VEC)):

    plt.plot(x,H_VEC[i,2]/sum(H_VEC[i,2]),'-',label='8, '+str(i))
    plt.xlabel('time in ns')

plt.legend()
plt.show()
"""
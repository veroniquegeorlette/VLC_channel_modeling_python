# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 10:47:22 2019

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

""" Paramètres de simulation """ 
C=3e8*1e-9 #vitesse de la lumière et en nanoseconde (m/ns)
rho= 0.8 #coefficient de reflexion
Ts=20#gain du filtre 
n=1.5#indice de réfraction



"""-----Emetteur----"""
 
theta = 50  #demi angle en degré 
m=-np.log10(2)/np.log10(np.cos(theta)) #Mode lambertien 
P_Total=1 #puissance totale transmise


"""----Recepteur----"""
Adet=np.exp(-4)#aire de detection du recepteur
FOV=60 #field of view
g=np.square(n)/np.square(np.sin(FOV))#gain 

"""----Dimensionnement de la piece ----"""
lx=5
ly=5
lz=3
#dimension de la pièce
Nx=lx*3
Ny=ly*3
Nz=np.round(lz*3)
#nombre de grille pour chaque surface
dA=(lz*ly)/(Ny*Nz)
#aire de la grille

x=np.linspace(int(-lx/2), int(lx/2), Nx)
y=np.linspace(int(-ly/2), int(ly/2), Ny)
z=np.linspace(-lz/2, lz/2, Nz)

xr, yr=np.meshgrid(x,y)

"""----Simulation--------"""
TP1=np.array([lx/4,ly/4,lz/2])
TP2=np.array([-lx/4,-ly/4,lz/2]) #Position de l'emetteur
TP3=np.array([lx/4,-ly/4,lz/2])
TP4=np.array([-lx/4,ly/4,lz/2])
H=(Nx, Ny)
h0=np.zeros(H)
h1=np.zeros(H)
h2=np.zeros(H)
h3=np.zeros(H)
h4=np.zeros(H)
delta_t=1/2
time_size=25
t_vector=np.arange(time_size/delta_t)
N_samp=int(time_size/delta_t)

H_VEC=np.zeros((Nx, Ny,len(t_vector)),float)

for ii in range(int(Nx)):
    
    for jj in range (int(Ny)):
        h_vector1= np.zeros(N_samp)
        h_vector2= np.zeros(N_samp)
        h_vector3= np.zeros(N_samp)
        h_vector4= np.zeros(N_samp)
        RP = np.array([x[ii], y[jj], -lz/2])#position de Rx
        """-------- TRANSMISSION LOS ----------""" 
        D1=np.sqrt(np.vdot(TP1-RP, TP1-RP))#Vecteur distance entre Tx et Rx
        tau1=D1/C
        index1=np.argwhere(round(tau1/delta_t)==t_vector)
        D2=np.sqrt(np.vdot(TP2-RP, TP2-RP))
        tau2=D2/C
        index2=np.argwhere(round(tau2/delta_t)==t_vector)
        D3=np.sqrt(np.vdot(TP3-RP, TP3-RP))
        tau3=D3/C
        index3=np.argwhere(round(tau3/delta_t)==t_vector)
        D4=np.sqrt(np.vdot(TP4-RP, TP4-RP))
        tau4=D4/C
        index4=np.argwhere(round(tau4/delta_t)==t_vector)
        cosphi=lz/D1
        cosphi2=lz/D2
        cosphi3=lz/D3
        cosphi4=lz/D4
        if np.abs(np.arccos(cosphi)*180/np.pi)<=FOV:
            TEST00=h0[ii, jj]
            TEST01=((m+1)*Adet*(cosphi**(m+1)))/(2*np.pi*np.square(D1))
            h0[ii, jj]= TEST00 + TEST01
            P_rec_A1=h0*P_Total*Ts*g
            h_vector1[index1]=h_vector1[index1]+((m+1)*Adet*(cosphi**(m+1)))/(2*np.pi*np.square(D1))*P_Total*Ts*g
            
        if np.abs(np.arccos(cosphi2)*180/np.pi)<=FOV:
            TEST00=h0[ii, jj]
            TEST02=((m+1)*Adet*(cosphi2**(m+1)))/(2*np.pi*np.square(D2))
            h0[ii, jj]= TEST00+TEST02
            P_rec_A1=h0*P_Total*Ts*g
            h_vector2[index2]=h_vector2[index2]+((m+1)*Adet*(cosphi2**(m+1)))/(2*np.pi*np.square(D2))*P_Total*Ts*g
            
        if np.abs(np.arccos(cosphi3)*180/np.pi)<=FOV:
            TEST00=h0[ii, jj]
            TEST03=((m+1)*Adet*(cosphi3**(m+1)))/(2*np.pi*np.square(D3))
            h0[ii, jj]= TEST00+TEST03
            P_rec_A1=h0*P_Total*Ts*g
            h_vector3[index3]=h_vector3[index3]+((m+1)*Adet*(cosphi3**(m+1)))/(2*np.pi*np.square(D3))*P_Total*Ts*g
           
        if np.abs(np.arccos(cosphi4)*180/np.pi)<=FOV:
            TEST00=h0[ii, jj]
            TEST04=((m+1)*Adet*(cosphi4**(m+1)))/(2*np.pi*np.square(D4))
            h0[ii, jj]= TEST00+TEST04
            P_rec_A1=h0*P_Total*Ts*g
            h_vector4[index4]=h_vector4[index4]+((m+1)*Adet*(cosphi4**(m+1)))/(2*np.pi*np.square(D4))*P_Total*Ts*g
        H_VEC[ii,jj]=np.asarray(h_vector1+h_vector2+h_vector3+h_vector4)
        

"""This is the classical Plot of the power received on the floor"""
H=h1+h2+h3+h4+h0
P_rec=H*P_Total*Ts*g

fig = plt.figure(figsize=plt.figaspect(2.))
ax = Axes3D(fig)
im=ax.plot_surface(xr, yr, P_rec, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
fig.colorbar(im,orientation="vertical")
plt.ylabel('y (m)')
plt.xlabel('x (m)')
plt.title('Received power')


fig.legend()
#print(H_VEC)

#%% 
"""This is thecomputation of the impulse response of the LOS contribution"""
fig = plt.figure(figsize=plt.figaspect(2.))
RX1=np.array([lx/4,ly/4])
H_0=P_rec[1,1]
x=np.arange(0,25,0.5)
""""
for i in range (14):
    for j in range(14):
        plt.figure()
        plt.plot(t_vector,H_VEC[i,j])
"""

#plt.plot(x,H_VEC[4,4],label='4,4')
#plt.plot(x,H_VEC[12,4],label='12,4')
##plt.plot(x,H_VEC[4,12],label='4,12')
#plt.plot(x,H_VEC[12,12],label='12,12')
#plt.plot(x,H_VEC[7,7]/sum(H_VEC[7,7]),label='Centre of the room')
plt.stem(x,H_VEC[3,4]/sum(H_VEC[3,4]),label='3,4, random spot to see the 4 contributions')
plt.ylabel('Normalized power contribution')
plt.xlabel('time in ns')
plt.title('Channel Impulse Response of a VLC system')
plt.legend(prop={'size': 8})
plt.show()
       

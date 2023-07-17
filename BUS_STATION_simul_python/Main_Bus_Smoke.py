# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 16:01:17 2019

@author: VeroGeorl
"""

"""
For the indoor function:
    inputs are : (FOV_, rho1_, rho2_, rho3_, rho4_, Ts_, n_, theta_, P_Total_, Adet_, lx_, ly_, lz_, Ng, TP1x, TP1y )
    outputs are: (P_rec, xr, yr, H, H_LOS, H_nLOS)
    
For the outdoor function:
    inputs are: (Vis,Lambda)
    output is: attenuation in dB/km
"""

#%%
"""Import section"""
import fINDOOR_Bus_Smoke as fin
import Plot_levelcurve_Bus_Smoke as plc
import fAttenuationsBus_Smoke as f_att
import numpy as np

#%%

""" Parameters section for the bus station under smoke or weather constraints"""
"""Emitter"""
#Half power angle(a 3dB)
theta=60
#Power emitted in mwatt
P=2000  
#Position of the emitter
Tx_pos=[0,0,1.54/2]
#Wavelength in nm
lam=450 #nm blue 
"""Receiver"""
#We suppose here the FDS100 Thorlabs
FOV=60 #Field of View
#Area=0.13#recepteur in cm²
Area=1.3E-05 # in m²
Ts=1 #Band pass filter
n=1.5#Refractive Index of the photoreceptor's lens
R=0.65 #A/W Responsivity
"""Room parameters"""
#dimension of the room
Space_=[2.9,2.1,1.54] 
#Reflexion coefficient of each wall. Some walls react differently for IR and VL
rho=[0.34 ,0.1,0.034,0] 
#Parameter linked to the number of divisions are made in the floor of the room
#Ng=20
Ng=10
 #Wall or not? Check image to see which wall corresponds to which index
Wall=[True,False,False,False]


#vis=np.logspace(np.log10(0.5),np.log10(60),50)#Visibility in km
vis=[0.3]


"""Other important parameters to compute the noise power"""
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





#%%
"""Computation of the power received in the receiver plane for the bus shelter"""


vi=0.3 #Visibility in km

P_rec, xr, yr, H, Hlos, Hnlos = fin.fINDOOR(FOV, rho, Ts, n, theta, P,Area, Space_, Ng, Tx_pos, Wall )
A_fog,A_fog_dB=f_att.compute_attenuation_fog(vi,lam)
A_smoke,A_smoke_dB=f_att.compute_attenuation_smoke_v(vi)
#P_r=10*np.log10(P_rec*H/1000)
P_r=10*np.log10(P_rec*H)
P_r_lin=P_rec*H
H_gaindB=10*np.log10(H)
            
plc.Save_data(xr, yr,P_r,theta,Wall,vi)
#plc.Plot_levelGain(xr, yr,H_gaindB,thet,w,vi,'level_gain')
#plc.Plot_level_lin(xr, yr,P_r_lin,thet,w,vi,'level_lin')



#%%
""" Computation of the noise """
"""Shot noise"""
# also related to background noise from paper on outdoor VLC for ITS
P_amb=bg_irra*BW_optF*Area*I2
Shot_noise=2*e*R*(P_r_lin*np.exp(-(A_fog+A_smoke)*Space_[2]*(10**(-3)))+P_amb)*D_r

"""
Thermal noise
"""
Term_1=(8*np.pi*kb*T*Capa*Area*I2*D_r**2)/G
Term_2=(16*np.pi**2*kb*T*GG*BW_optF**2*Area**2*I3*D_r**3)/gm
Thermal_noise=Term_1+Term_2


#%%
"""
Plot the data saved in the .dat files 
"""
#x=np.linspace(-Space_[0]/2, Space_[0]/2, int (Space_[0]*Ng))
#y=np.linspace(-Space_[1]/2, Space_[1]/2, int (Space_[1]*Ng))
x=np.linspace(0, Space_[0], int (Space_[0]*Ng))
y=np.linspace(0, Space_[1], int (Space_[1]*Ng))
xr, yr=np.meshgrid(x,y)


data=[]
data.append(np.loadtxt( r'P_r_60_'+str(Wall)+'_'+str(round(vi,3))+'.dat' ))
        
Min_all=np.min(data)
Max_all=np.max(data)
i=0
plc.Just_plot(xr,yr,Min_all,Max_all,data[i].transpose(),rho,Wall,i)


#%%
"""
Computation of the SNR
"""
""" 
This part computes the SNR of the transmission in visible light
"""
#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

SNR=[]
Num=(R*P_r_lin)**2.0
Num_att=(R*P_r_lin*np.exp(-(A_fog+A_smoke)*Space_[2]*(10**(-3))))**2.0
Den=Shot_noise+Thermal_noise
SNR=Num/Den
SNR_=Num_att/Den
"""
fig = plt.figure(figsize=plt.figaspect(2.))
ax = Axes3D(fig)
surf = ax.plot_surface(xr, yr, (10*np.log10(SNR)).transpose(), rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
#surf = ax.plot_surface(xr, yr, SNR.transpose(), rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
#ax.set_zlim(v_min, v_max)
#ax.zaxis.set_major_locator(LinearLocator(5))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('SNR')
plt.title("Signal to Noise Ratio dB ")
fig.colorbar(surf, shrink=0.3, aspect=10)
plt.show()
"""

fig = plt.figure(figsize=plt.figaspect(2.))
ax = Axes3D(fig)
surfi = ax.plot_surface(xr, yr, (10*np.log10(SNR-SNR_)).transpose(), rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
#surfi = ax.plot_surface(xr, yr, (10*np.log10(SNR-SNR_)).transpose(), rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
#surf = ax.plot_surface(xr, yr, SNR.transpose(), rstride=1, cstride=1, cmap='jet',linewidth=0, antialiased=False)
#ax.set_zlim(v_min, v_max)
#ax.zaxis.set_major_locator(LinearLocator(5))
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Δ SNR (dB)')
plt.title("SNR difference when there is or not attenuation due to smoke and fog ")
fig.colorbar(surfi, shrink=0.3, aspect=10)
plt.show()
 


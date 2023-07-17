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
import fINDOOR_House_Smoke as fin
import Plot_levelcurve_House_Smoke as plc
import fAttenuationsHouse_Smoke as f_att
import numpy as np

#%%

""" Parameters section for the bus station under smoke or weather constraints"""
"""Emitter"""
#Half power angle(a 3dB)
theta=60
#Power emitted in mwatt
#P=27000 
P=27 #in watt 
#Position of the emitter
Tx_pos=[-2.5,0,3]



"""Receiver"""
FOV=60 #Field of View
#Area=0.13#recepteur in cm²
Area=1.3E-05 #recepteur in m²
Ts=1 #Band pass filter
n=1.5#Refractive Index of the photoreceptor's lens


"""Room parameters"""
#dimension of the room
Space_=[6,6,6] 
#Reflexion coefficient of each wall. Some walls react differently for IR and VL
rho=[0.34 ,0,0,0] 
#Parameter linked to the number of divisions are made in the floor of the room
Ng=3
#Wall or not? Check image to see which wall corresponds to which index
Wall=[True,False,False,False]
#Wavelength in um
lam=0.650 



#vis=np.logspace(np.log10(0.5),np.log10(60),50)#Visibility in km
vis=[0.3]








#%%
"""Tests for the ideal case with the position of the emitter changing"""


vi=0

P_rec, xr, yr, H, Hlos, Hnlos = fin.fINDOOR(FOV, rho, Ts, n, theta, P,Area, Space_, Ng, Tx_pos, Wall )

#P_r=10*np.log10(P_rec*H/1000)
P_r=10*np.log10(P_rec*H)
P_r_lin=P_rec*Hlos
H_gaindB=10*np.log10(H)
            
plc.Save_data(xr, yr,P_r,theta,Wall,vi)
#plc.Plot_levelGain(xr, yr,H_gaindB,thet,w,vi,'level_gain')
#plc.Plot_level_lin(xr, yr,P_r_lin,thet,w,vi,'level_lin')


#%%
"""
Plot the data saved in the .dat files 
"""
x=np.linspace(0, Space_[0], int (Space_[0]*Ng))
y=np.linspace(0, Space_[1], int (Space_[1]*Ng))
xr, yr=np.meshgrid(x,y)


data=[]
data.append(np.loadtxt( r'P_r_60_'+str(Wall)+'_'+str(round(vi,3))+'.dat' ))
        
Min_all=np.min(data)
Max_all=np.max(data)
i=0
plc.Just_plot(xr,yr,Min_all,Max_all,data[i].transpose(),rho,Wall,i)


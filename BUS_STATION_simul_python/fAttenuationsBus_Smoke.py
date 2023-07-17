# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:49:15 2019

@author: VeroGeorl
Equations taken from "Analysis of Fog Effects on Terrestrial FSO communication Links
First the Kim model was used and then the Ijaz model for the smaller visibility
"""
import numpy as np

#For the parameters, Vis should be in kilometers and lambda in micrometers

def compute_attenuation_fog(Vis,Lamb):
    att=0.0
    a=0.0
    Vis_=int(round(Vis))
    #Vis_='%.3f'%(Vis)
    if(0.600<Lamb<1.550): #wavelength in um
        if(Vis>=50):
            q = 1.6
            a = ((3.91/Vis)*((Lamb/0.550)**(-q)))
            att = 10*a/np.log(10)
            #print('1')
        elif(6<=Vis<50):
            q = 1.3
            a = ((3.91/Vis)*((Lamb/0.550)**(-q)))
            att = 10*a/np.log(10)
            #print('2')
        elif(1<=Vis<6):
            q = (0.16*Vis)+0.34
            a = ((3.91/Vis)*((Lamb/0.550)**(-q)))
            att = 10*a/np.log(10)
            #print('3')
        elif(0.5<Vis<1):
             q = Vis-0.5
             a = ((3.91/Vis)*((Lamb/0.550)**(-q)))
             att = 10*a/np.log(10)
             #print('4')
        elif(0.015<Vis<0.5): #Ijaz model
            q=0.1428*Lamb-0.0947 #Fog equation
            a=((17/Vis)*((Lamb/0.550)**(-q)))
            att = 10*a/np.log(10)
            #print('5')
        elif(Vis_<0.15):
            print("Attenuation approximated to 315 dB/km")
            a=72.53
            att=315
            #print('6')
    else: print("Wavelength out of range")
    return att, a
"""
def compute_attenuation_smoke(Vis,Lamb):
    att=0.0
    q=0.8467*Lamb-0.5212 #Smoke equation
    a=(17/Vis)*(Lamb/550)**(-q)
    att = 10*a/np.log(10)
    return att
"""
def compute_attenuation_smoke(Lamb,dimensions,mass):
    sigma=((632.8/Lamb)*8500)/1000 #m²/g or 7.6
    #Km=np.random.normal(8.700,1.100,1)
    Volume=dimensions[0]*dimensions[1]*dimensions[2]
    epsilon={"Douglas fir"   :  np.random.normal(0.0175,0.0075,1),
             "Hardboard"    : 	 np.random.normal(0.0007,0.0003,1),
             "Fiberboard"   :  	np.random.normal(0.0075,0.0025,1),
             "Polyvinylchloride "  :  	0.12,
             "Polyurethane (flexible)":	np.random.normal(0.0225,0.0125,1),
             "Polyurethane (rigid)" : 	0.09,
             "Polystyrene "   :	np.random.normal(0.15,0.02,1) ,
             "Polypropylene " : 	np.random.normal(0.008,0.008,1),
             "Polymethylmethacrylate"  :	0.02 ,
             "Polyoxymethylene" :	0,
             "Cellulosic" :	0.015 }
    eps=list(epsilon.values())
    eps_mean=np.mean(eps)
    m=mass
    att_lin=sigma*m*eps_mean/Volume
    att_dB_km = 10*att_lin/np.log(10)
    return att_lin, att_dB_km

def compute_attenuation_smoke_v(Visibility):
    Att_coef_s = 0.0
    Att_coef_s_dB_km=0.0
    Att_coef_s=3/(Visibility)
    Att_coef_s_dB_km=10*Att_coef_s/np.log(10)
    return Att_coef_s,Att_coef_s_dB_km
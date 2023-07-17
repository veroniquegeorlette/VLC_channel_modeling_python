# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:28:20 2019
Updated 02/06/2020 15:00

@author: VeroGeorl
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import numpy as np
import matplotlib.pyplot as plt
import os

from mpl_toolkits.mplot3d import Axes3D

"""
def Just_plot(xr,yr,P_r,thet,wa,vi) :
    plt.figure(figsize=(16,6))
    plt.subplot(1,2,1,projection='3d')
    plt.gca(projection='3d').plot_surface(xr, yr,P_r,rstride=1,cstride=1,linewidth=0,
    cmap='jet',antialiased=False)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Puissance reçue (dBm)')
    plt.subplot(1,2,2)
    vec=np.linspace(np.min(P_r),np.max(P_r),10)
    cont=plt.contour(xr,yr,P_r,vec,cmap='jet')
    plt.clabel(cont)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Courbes de niveau Vis: '+str(vi)) 
    plt.show()
    return None
"""
def Just_plot(xr,yr,min_,max_,P_r,thet,wa,vi) :
    plt.figure(figsize=(16,6))
    plt.subplot(1,2,1,projection='3d')
    ax=plt.gca(projection='3d')
    ax.plot_surface(xr, yr,P_r,rstride=1,cstride=1,vmin=min_,vmax=max_,linewidth=0,
    cmap='jet',antialiased=False)
    ax.set_zlim(min_, max_)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Power received (dBm)')
    #plt.axis('equal')
    #plt.axis("equal")
    plt.subplot(1,2,2)
    vec=np.linspace(min_,max_,10)
    cont=plt.contour(xr,yr,P_r,vec,cmap='jet')
    plt.clabel(cont)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Level curves dBm') 
    
    #plt.axis("equal")
    #plt.axis('equal')
    plt.rc('font', size=14) #controls default text size
    plt.rc('axes', titlesize=14) #fontsize of the title
    plt.rc('axes', labelsize=14) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=14) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=14) #fontsize of the y tick labels
    plt.rc('legend', fontsize=14) #fontsize of the legend
    plt.show()
    plt.draw()
    return None

def Plot_levelcurves(xr,yr,P_r,thet,wa,vi,dirName) :
    plt.figure(figsize=(16,6))
    plt.subplot(1,2,1,projection='3d')
    plt.gca(projection='3d').plot_surface(xr, yr,P_r,rstride=1,cstride=1,linewidth=0,
    cmap='jet',antialiased=False)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Puissance reçue (dB)')
    plt.subplot(1,2,2)
    vec=np.linspace(np.min(P_r),np.max(P_r),10)
    cont=plt.contour(xr,yr,P_r,vec,cmap='jet')
    plt.clabel(cont)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Courbes de niveau '+str(vi)) 
    Vi_='%.3f'%(vi)
    
 
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
        
    path=dirName+r'\fi'+str(vi)+'.png'
    plt.savefig(path)
    plt.close()
    return None

def Plot_level_lin(xr,yr,P_r,thet,wa,vi,dirName) :
    plt.figure(figsize=(16,6))
    plt.subplot(1,2,1,projection='3d')
    plt.gca(projection='3d').plot_surface(xr, yr,P_r,rstride=1,cstride=1,linewidth=0,
    cmap='jet',antialiased=False)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Puissance reçue (W)')
    plt.subplot(1,2,2)
    vec=np.linspace(np.min(P_r),np.max(P_r),10)
    cont=plt.contour(xr,yr,P_r,vec,cmap='jet')
    plt.clabel(cont)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Courbes de niveau '+str(vi)) 
    Vi_='%.3f'%(vi)
    
 
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
 
   
    path=dirName+r'\fi'+str(vi)+'.png'
    plt.savefig(path)
    plt.close()
    return None

def Plot_levelGain(xr,yr,P_r,thet,wa,vi,dirName) :
    plt.figure(figsize=(16,6))
    plt.subplot(1,2,1,projection='3d')
    plt.gca(projection='3d').plot_surface(xr, yr,P_r,rstride=1,cstride=1,linewidth=0,
    cmap='jet',antialiased=False)
    plt.xlabel('x (m)',fontsize=14)
    plt.ylabel('y (m)')
    plt.title('Channel DC Gain  dB')
    plt.subplot(1,2,2)
    vec=np.linspace(np.min(P_r),np.max(P_r),10)
    cont=plt.contour(xr,yr,P_r,vec,cmap='jet')
    plt.clabel(cont)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Level curve channel gain '+str(vi)) 
    Vi_='%.3f'%(vi)
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
 
   
    path=dirName+r'\fi'+str(vi)+'.png'
    plt.savefig(path)
    plt.close()
    return None


def Save_data(xr,yr,P_r,thet,wa,vi):
    np.savetxt('xr_'+'_'+str(thet)+'_'+str(wa)+'_'+str(vi)+'.dat', xr)
    np.savetxt('yr_'+'_'+str(thet)+'_'+str(wa)+'_'+str(vi)+'.dat', yr)
    np.savetxt('P_r'+'_'+str(thet)+'_'+str(wa)+'_'+str(vi)+'.dat', P_r)
    return None

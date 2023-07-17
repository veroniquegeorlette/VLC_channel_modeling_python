# -*- coding: utf-8 -*-
"""
Created on Fri May  7 21:34:49 2021

@author: VeroGeorl
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def get_cartesian(theta,phi,rho):
        #precalculate sine and cosine of phi and theta
        cos_phi = np.cos(phi);
        sin_phi = np.sin(phi);
        cos_theta = np.cos(np.pi-theta);
        sin_theta = np.sin(np.pi-theta);
        
        x = rho *sin_theta *cos_phi
        y = rho *sin_theta *sin_phi
        z = rho *(cos_theta)
        #reshape to flat (n x 1)vectors:
        x = x.flatten()
        y = y.flatten()
        z = z.flatten()
        #end by returning point coordinates
        return x,y,z

with open('test_indus.txt','r') as file:
    indus = file.read()


ind=indus.strip()


hey=" ".join(ind.split())


final=hey.replace(" ", ",")

pattern=list(map(float,final.split(',')))


rho_=np.asarray(pattern)
rho=rho_.reshape(19,42)


theta=np.arange(0,102.5,2.5) 
np.append(theta,180)
theta=theta* np.pi / 180.0
phi=np.arange(0,95,5)* np.pi / 180.0

coord=[]
for i in range(phi.size-1):
    for j in range(theta.size-1):
        #print(rho[i][j])
        x,y,z=get_cartesian(theta[j],phi[i],rho[i][j])
        coord.append([x,y,z])
        
A=np.asarray(coord)
B=A.transpose()
x_vec=B[0][0]
y_vec=B[0][1]
z_vec=B[0][2]
X,Y = np.meshgrid(x_vec,y_vec)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


pattern = ax.scatter(x_vec,y_vec, z_vec,c=z_vec, linewidth=0, antialiased=False)
#ax.plot([0,0], [0,0],zs=[0,np.max(z_vec)])

fig.colorbar(pattern,orientation="vertical")
plt.show()

#%%
"""
This part of the code will compute the position of the receiver with respect
to the emitter in spherical coordinates. 
The center of the axes is located on the light.

"""
efficacite_lum=130 #lum/sr
Area_rx=1.30e-05 #m²

def get_spherical(x,y,z):
    rho=np.sqrt(x**2+y**2+z**2)
    phi=np.arctan2(y,x)
    if phi<0:
        phi=phi+2*np.pi
    #the minus on the z is to have the same reference angle as the IES file
    theta=np.arccos(z/np.sqrt(x**2+y**2+z**2))
    return rho,phi*180/np.pi,180-theta*180/np.pi

TX_pos=np.array([0,0,0])
RX_pos=np.array([1,1,-2])

RX_pos_spher=get_spherical(RX_pos[0],RX_pos[1],RX_pos[2])
print("RX_pos_spher",RX_pos_spher)

for_theta=int(round(RX_pos_spher[2]/2.5))
for_phi=int(round(RX_pos_spher[1]/5))
Puissance_lu=rho[for_phi][for_theta]

print(Puissance_lu)
Puissance_W_sr=Puissance_lu/efficacite_lum
print(Puissance_W_sr)

P_rx_W=Puissance_W_sr*Area_rx/(RX_pos_spher[0]**2)
print(P_rx_W)
P_rx_dBm=10*np.log10(P_rx_W/0.001)
print(P_rx_dBm)



#%%
"""
This part will build a plane in 3D within the space defined above to assess
the power reaching the receiver when it is located at 1.5 m above the ground.
An interpolation function is used to smooth the data on the phi angle 
"""
def get_spherical(x,y,z):
    rho=np.sqrt(x**2+y**2+z**2)
    phi=np.arctan2(y,x)
    if phi<0:
        phi=phi+2*np.pi
    #the minus on the z is to have the same reference angle as the IES file
    theta=np.arccos(z/np.sqrt(x**2+y**2+z**2))
    return rho,phi*180/np.pi,180-(theta*180/np.pi)

#This function returns the power per solid angle from the table of 
#the radiation pattern
def get_W_per_sr(RX_pos):
    RX_pos_spher=get_spherical(RX_pos[0],RX_pos[1],RX_pos[2])
    for_theta=int(round(RX_pos_spher[2]/2.5))
    Puissance_lu=np.interp((RX_pos_spher[1])/5,phi[:-1]*180/np.pi,Rho_[:,for_theta])
    Puissance_W_sr=Puissance_lu/efficacite_lum
    return Puissance_W_sr

efficacite_lum=130 #lum/sr
Area_rx=1.30e-05 #m²
TX_pos=np.array([0,0,0])
#RX_pos=np.array([1,1,-2])
n=1.5
FOV=60
FOV_rad=FOV*np.pi/180
g=np.square(n)/np.square(np.sin(FOV_rad))#gain 

z_pos=-1.5

H=(len(x),len(y))
Rho_=np.delete(rho, -1, 0)
#h=np.zeros(H)
h=[]
h_dBm=[]

lx=3
ly=3
lz=3 # m hauteur de la lampe
Nx=50
Ny=50
Nz=100 #nombre de case à évaluer sur mon mur
x=np.linspace(0,lx,Nx)
y=np.linspace(0,ly,Ny)
z=np.linspace(-lz/2,lz/2,Nz)
for i in x:
    for j in y:
        RX_pos=np.array([i, j, z_pos])
        
        
        RX_pos_spher=get_spherical(RX_pos[0],RX_pos[1],RX_pos[2])
        #here=180-RX_pos_spher[2]
        
        if(RX_pos_spher[2]<=FOV):
            #for_phi2=int(round(RX_pos_spher[1])/30.)
            #for_theta=int(round(RX_pos_spher[2]))
            #Puissance_lu=np.interp((RX_pos_spher[1])/30.,phi[:-1]*180/np.pi,Rho_[:,for_theta])
            #print( 'phi  '+str(for_phi)+ 'theta  ' + str(for_theta))
            #Puissance_lu1=rho[for_phi][for_theta]
            #Puissance_W_sr=Puissance_lu/efficacite_lum
            
            Puissance_W_sr=get_W_per_sr(RX_pos)
            P_rx_W=Puissance_W_sr*Area_rx/(RX_pos_spher[0]**2)
            P_rx_W_g=P_rx_W*g
            h.append(P_rx_W_g)
            P_rx_dBm=10*np.log10(P_rx_W_g/0.001)
            h_dBm.append(P_rx_dBm)
            #TEST00=H[ii, jj]
            #TEST01=((m+1)*Adet*(cosphi**(m)))/(2*np.pi*np.square(D1))

        else:
            h.append(0)
            h_dBm.append(np.log10(0/0.001))

       
X,Y=np.meshgrid(x,y)
#Power in Watts
h_array=np.asarray(h)  
h_matrix=h_array.reshape(len(x),len(y))
#Power in dBm
h_array_dBm=np.asarray(h_dBm)   
h_matrix_dBm=h_array_dBm.reshape(len(x),len(y))


ax = fig.add_subplot(111, projection='3d')
plt.figure(figsize=(16,6))
plt.subplot(1,2,1,projection='3d')
ax=plt.gca(projection='3d')
ax.plot_surface(X,Y, h_matrix*1000, cmap='jet')
plt.title('Power reaching the receiver plane in mW') 
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.subplot(1,2,2)
vec=np.linspace(np.min(h_matrix*1000), np.max(h_matrix*1000),10)
cont=plt.contour(X, Y, h_matrix*1000,vec,cmap='jet')
plt.clabel(cont)
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Level curves mW') 
#plt.axis('equal')
plt.show()



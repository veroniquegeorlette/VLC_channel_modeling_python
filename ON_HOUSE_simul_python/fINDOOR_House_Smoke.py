# -*- coding: utf-8 -*-
import numpy as np

def fINDOOR(FOV_, rho, Ts_, n_, theta_, P_Total_, Adet_, Space_, Ng, Tx_pos,wall ):
    
    """ Paramètres de simulation """ 

    rho1= rho[0] #coefficient de reflexion de chaque mur 
    rho2= rho[1]
    rho3= rho[2]
    rho4= rho[3]
    wall_1=wall[0]
    wall_2=wall[1]
    wall_3=wall[2]
    wall_4=wall[3]
    
    Ts=Ts_#gain du filtre    
    
    n=n_#indice de réfraction

    #c=300000

    """-----Emetteur----"""
 
    theta = theta_*(np.pi/180)  #demi angle en degré 
    m=-np.log10(2)/np.log10(np.cos(theta)) #Mode lambertien 
    P_Total=P_Total_ #puissance totale transmise


    """----Recepteur----"""

    Adet=Adet_#aire de detection du recepteur
    FOV=FOV_*(np.pi/180) #field of view
    FOV_deg=FOV_
    g=np.square(n)/np.square(np.sin(FOV))#gain 

    """----Dimensionnement de la piece ----"""
    
    lx= Space_[0]
    ly= Space_[1]
    lz= Space_[2]
    #dimension de la pièce
    Nx=int (lx*Ng)
    Ny=int (ly*Ng)
    Nz=np.round(lz*Ng)
    #nombre de grille pour chaque surface
    dA=(lz*ly)/(Ny*Nz)
    #aire de la grille
    
    x=np.linspace(-lx/2, lx/2, Nx)
    y=np.linspace(-ly/2, ly/2, Ny)
    z=np.linspace(-lz/2, lz/2, Nz)

    xr, yr=np.meshgrid(x,y)

    """----Simulation--------"""
    
    TP1=np.array([Tx_pos[0],Tx_pos[1],Tx_pos[2]]) #Position de l'emetteur
    H=(Nx, Ny)

    
    H0=np.zeros(H)
   # h0=np.zeros(H)
    
    H1=np.zeros(H)
   # h1=np.zeros(H)
    
    H2=np.zeros(H)
   # h2=np.zeros(H)
    
    H3=np.zeros(H)
   # h3=np.zeros(H)
    
    H4=np.zeros(H)
   # h4=np.zeros(H)
    
    for ii in range(int(Nx)):
        
        for jj in range (int(Ny)):
        
            RP = np.array([x[ii], y[jj], -lz/2])#position de Rx
        
            """-------- TRANSMISSION LOS ----------""" 
        
            D1=np.sqrt(np.vdot(TP1-RP, TP1-RP))#Vecteur distance entre Tx et Rx
            cosphi=lz/D1
            
            if np.abs(np.arccos(cosphi))<=FOV:
                TEST00=H0[ii, jj]
                TEST01=((m+1)*Adet*(cosphi**(m)))/(2*np.pi*np.square(D1))
                H0[ii, jj]= TEST00 + TEST01
                
                
            
            """------- REFLEXION MUR 1 -----------"""
            if(wall_1==True):
                for kk in range (int(Ny)):
                
                    for ll in range (int(Nz)):
                        WP1=np.array([-lx/2, y[kk], z[ll]]) #point de l'incidence
                        D1=np.sqrt(np.vdot(TP1-WP1, TP1-WP1))
                    # distance entre Tx et point d'incidence
                        cos_phi=(np.abs(WP1[2]-TP1[2]))/D1
                        cos_alpha=(np.abs(TP1[0]-WP1[0]))/D1
                    
                        D2=np.sqrt(np.vdot(WP1-RP, WP1-RP))
                    #distance entre le point d'incidence et Rx
                        cos_psy=(np.abs(WP1[2]-RP[2]))/D2
                        cos_beta=(np.abs(WP1[0]-RP[0]))/D2
                
                        if np.abs(np.arccos(cos_psy))<=FOV:
    
                        
                            TEST10= H1[ii, jj]
                            TEST11= ((m+1)*Adet*rho1*dA*(cos_phi**(m))
                            *cos_alpha*cos_beta*cos_psy)/(2*np.square(np.pi)*np.square(D1)*np.square(D2))
                        
                            H1[ii, jj]= TEST10 + TEST11
                        
                  
            """------- REFLEXION MUR 2 ----------"""
            if(wall_2==True):
                for kk in range (int(Nx)):
                
                    for ll in range (int(Nz)):
                        WP2=np.array([x[kk], -ly/2, z[ll]]) #point de l'incidence
                        D1=np.sqrt(np.vdot(TP1-WP2, TP1-WP2))
                    # distance entre Tx et point d'incidence
                        cos_phi=(np.abs(WP2[2]-TP1[2]))/D1
                        cos_alpha=(np.abs(TP1[1]-WP2[1]))/D1
                    
                        D2=np.sqrt(np.vdot(WP2-RP, WP2-RP))
                    #distance entre le point d'incidence et Rx
                        cos_psy=(np.abs(WP2[2]-RP[2]))/D2
                        cos_beta=(np.abs(WP2[1]-RP[1]))/D2
                    
                        if np.abs(np.arccos(cos_psy))<=FOV:
    
                            TEST20= H2[ii, jj]
                            TEST21= ((m+1)*Adet*rho2*dA*(cos_phi**(m))
                            *cos_alpha*cos_beta*cos_psy)/(2*np.square(np.pi)*np.square(D1)*np.square(D2))
                            H2[ii, jj]= TEST20 + TEST21
                        
            """------- REFLEXION MUR 3  ----------"""
            if(wall_3==True):
                for kk in range (int(Ny)):
                
                    for ll in range (int(Nz)):
                        WP3=np.array([lx/2, y[kk], z[ll]]) #point de l'incidence
                        D1=np.sqrt(np.vdot(TP1-WP3, TP1-WP3))
                    # distance entre Tx et point d'incidence
                        cos_phi=(np.abs(WP3[2]-TP1[2]))/D1
                        cos_alpha=(np.abs(TP1[0]-WP3[0]))/D1
                    
                        D2=np.sqrt(np.vdot(WP3-RP, WP3-RP))
                    #distance entre le point d'incidence et Rx
                        cos_psy=(np.abs(WP3[2]-RP[2]))/D2
                        cos_beta=(np.abs(WP3[0]-RP[0]))/D2
                    
                        if np.abs(np.arccos(cos_psy))<=FOV:
    
                        
                            TEST30= H3[ii, jj]
                            TEST31= ((m+1)*Adet*rho3*dA*(cos_phi**(m))
                            *cos_alpha*cos_beta*cos_psy)/(2*np.square(np.pi)*np.square(D1)*np.square(D2))
                        
                            H3[ii, jj]= TEST30 + TEST31
                   
            """------- REFLEXION MUR 4 ----------"""
            if(wall_4==True):
                for kk in range (int(Nx)):
                
                    for ll in range (int(Nz)):
                        WP4=np.array([x[kk], ly/2, z[ll]]) #point de l'incidence
                        D1=np.sqrt(np.vdot(TP1-WP4, TP1-WP4))
                    # distance entre Tx et point d'incidence
                        cos_phi=(np.abs(WP4[2]-TP1[2]))/D1
                        cos_alpha=(np.abs(TP1[1]-WP4[1]))/D1
                    
                        D2=np.sqrt(np.vdot(WP4-RP, WP4-RP))
                    #distance entre le point d'incidence et Rx
                        cos_psy=(np.abs(WP4[2]-RP[2]))/D2
                        cos_beta=(np.abs(WP4[1]-RP[1]))/D2
                    
                        if np.abs(np.arccos(cos_psy))<=FOV:
    
                            TEST40= H4[ii, jj]
                            TEST41= ((m+1)*Adet*rho4*dA*(cos_phi**(m))
                            *cos_alpha*cos_beta*cos_psy)/(2*np.square(np.pi)*np.square(D1)*np.square(D2))
                            H4[ii, jj]= TEST40 + TEST41
    
    
    H_LOS=H0
    H_nLOS=H1+H2+H3+H4
    H= H_LOS+H_nLOS
    P_rec=P_Total*Ts*g
    #P_rec_dB = 10*np.log10(P_rec)
  

    
    return(P_rec, xr, yr, H, H_LOS, H_nLOS)



        





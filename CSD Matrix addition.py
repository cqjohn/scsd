#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 22:07:08 2023

@author: genah
"""

#import libraries 
import matplotlib.pyplot as plt
import numpy as np 



#define constants 
e = 1 #normalised everything so not 1.6e-19
pd = 4000 #pd 

#charging energies and coupling energy 

def Ec1func(C1,C2,Cm):
    cmsq = Cm **2
    frac = cmsq/(C1*C2)
    val = ((e**2)/C1) * (1/(1-(frac)))
    return val 

def Ec2func(C1,C2,Cm):
    cmsq = Cm **2
    frac = cmsq/(C1*C2)
    val = ((e**2)/C2) * (1/(1-(frac)))
    return val 
    
def Ecmfunc(C1,C2,Cm):
    cmsq = Cm **2
    frac = (C1*C2)/cmsq
    val = ((e**2)/Cm) * (1/((frac)-1))
    return val 


#electrochemical potentials 

def mu1(N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2):
    val = (N1 - 0.5)*Ec1 + N2*Ecm - (1/e)*(Cg1*Vg1*Ec1 + Cg2*Vg2*Ecm)
    return val

def mu2(N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2):
    val = (N2 - 0.5)*Ec2 + N1*Ecm - (1/e)*(Cg1*Vg1*Ecm + Cg2*Vg2*Ec2)
    return val

#plot function 
def plotChargeStab(CL,CR,Cm,Cg1,Cg2):
    
    C1 = CL + Cm + Cg1
    C2 = CR + Cm + Cg2
    
    
    Vg1 = np.linspace(0,5,pd)#*1e-6
    Vg2 = np.linspace(0,5,pd)#*1e-6 
    
    Ec1 = Ec1func(C1,C2,Cm)
    Ec2 = Ec2func(C1,C2,Cm)
    Ecm = Ecmfunc(C1,C2,Cm)
    
    X, Y = np.meshgrid(Vg1, Vg2)
    fig, ax = plt.subplots()
    
    #make seperate array  for triple points and degenerate lines 
    #so that I can multiply triple point array by 2 to make it visible in plot
    Z = np.bool8(np.zeros([pd,pd]))
    Z_tri = np.bool8(np.zeros([pd,pd]))
    
    k = 1
    for i in range(0,4):
        for j in range(0,4):
            print(i,j)
            #STEP 1: Calculate mu1 and mu2 for i,j occupancy and adjacnet occupancy
            
            #calculate mu1 and mu2 for i,j where i and j are different electron occupancy 
            mu1_arr = np.around( np.array(mu1(i,j,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)), decimals=3) 
            mu2_arr = np.around( np.array(mu2(i,j,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)) , decimals=3) 
            
            #calculate adjacent mu1 and mu2 for electron transfer 
            mu1_arr_adj_i =  np.around( np.array(mu1(i+1,j,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)) , decimals=3) 
            mu2_arr_adj_j =  np.around( np.array(mu2(i,j+1,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)) , decimals=3) 
            
            
            mu1_arr_adj_j =  np.around( np.array(mu1(i,j+1,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)) , decimals=3) 
            mu2_arr_adj_i =  np.around( np.array(mu2(i+1,j,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)) , decimals=3) 
            
            #calculate opposite mu1 and mu2 for electron transfer 
            mu1_arr_opp =  np.around( np.array(mu1(i+1,j+1,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)) , decimals=3) 
            mu2_arr_opp =  np.around( np.array(mu2(i+1,j+1,Ec1,Ec2,Ecm,Cg1,Cg2,X,Y)) , decimals=3) 
            
            #store old values so not lost when we redefine them 
            Z_old = Z
            Z_old_tri = Z_tri
        
            #CONDITIONS ARE HERE -------------------------------------------->
        
            #triple point conditions 
            Z_new_tri_el = ((mu1_arr_adj_i == mu2_arr_adj_j) & (mu1_arr_adj_i ==0))
            Z_new_tri_hl = ((mu1_arr_opp == mu2_arr_opp) & (mu1_arr_opp ==0))
            
            
            #find the x and y coordiantes of the triple points 
            if np.where(Z_new_tri_el == 1)[0].size == 0:
                Z_1 = np.bool8(np.zeros([pd,pd]))
                Z_2 = np.bool8(np.zeros([pd,pd]))
                Z_3 = np.bool8(np.zeros([pd,pd]))
                Z_tri = Z_old_tri
            else:
                el_1, el_2 = np.where(Z_new_tri_el ==1)[0][0], np.where(Z_new_tri_el ==1)[1][0]
                if np.where(Z_new_tri_hl == 1)[0].size ==0:
                    Z_tri = Z_new_tri_el | Z_old_tri
                else:
                    Z_tri = Z_new_tri_el | Z_old_tri| Z_new_tri_hl
                    hl_1, hl_2 = np.where(Z_new_tri_hl ==1)[0][0], np.where(Z_new_tri_hl ==1)[1][0]
            
            
                    #Set Z_1 condition and restrict using triple point coordiantes 
                    Z_1 = (mu1_arr == mu2_arr)
                    Z_1[0:el_1,:] , Z_1[hl_1:,:], Z_1[:,0:el_2], Z_1[:,hl_2:]   = False , False, False, False
            
            
                #Find adjacent (hole) triple points for Z_2 and Z_3 conditions 
                hl_left =  ((mu1_arr_adj_j == mu2_arr_adj_j) & (mu1_arr_adj_j == 0))
                hl_right = ((mu1_arr_adj_i == mu2_arr_adj_i ) & (mu1_arr_adj_i == 0))
                

                
        
                
                #find the x and y coordiantes of the adjacent hole triple points 
                if np.where(hl_left == 1)[0].size == 0:
                    Z_2 = np.bool8(np.zeros([pd,pd]))
                    print("Adjacent left hole not found ")
                else: 
                    print("Ajacent left hole found")
                    hl_l_1, hl_l_2 = np.where(hl_left == 1)[0][0], np.where(hl_left == 1)[1][0]
                
                    #Set Z_2 condition and restrict using triple point coordiantes 
                    #left hole 
                    Z_2 = (mu2_arr_adj_j == 0)
                    #VG1 restrictions 
                    Z_2[:, 0:hl_l_2] = False 
                    Z_2[:, el_2:] = False 
                    #VG2 restrictions 
                    #Z_2[0:el_1, :] = False
                    #Z_2[hl_l_1:, 0] = False 
                    
                
                if np.where(hl_right == 1)[0].size == 0: 
                    Z_3 = np.bool8(np.zeros([pd,pd]))
                    print("Adjacent right hole not found ")
                else: 
                    print("Ajacent right hole found")
                    hl_r_1, hl_r_2 = np.where(hl_right == 1)[0][0], np.where(hl_right == 1)[0][1]
                    #Set Z_3 condition and restrict using triple point coordiantes 
                    #right hole 
                    Z_3 = (mu1_arr_adj_i == 0)
                    #VG1 restrictions 
                    #Z_3[:,0: el_2] = False  
                    #Z_3[:, hl_r_2:] = False 
                    #VG2 restrictions 
                    Z_3[0:hl_r_1,:] = False 
                    Z_3[el_1:,:] = False 
                            
            #CONDITIONS END HERE -------------------------------------------->
            
            #add three seperate conditions for degenerate lines
            Z_new =  Z_1 | Z_2 | Z_3
     
            #add old Z and new Z plots to update figure 
            Z =  Z_old  | Z_new
        
            #when plotting mutiply triple points to make distinct from lines 
            Z_plot =  (Z_tri*2 + Z + Z_new_tri_el)*4
            
            
            #standard plot code 
            contourplot = plt.contourf(X,Y,Z_plot, cmap=plt.cm.bone,
                  origin='lower')
            ax.set_aspect('equal')
            cbar = plt.colorbar(contourplot)
            plt.show()
    np.save("data.npy",Z_plot )
            

        
        

    
#CL, CR, Cm, Cg1, Cg2   
#change third parameter to Cm = 0.00000000001 for uncoupled plot   
plotChargeStab(0.6,0.5,0.2,0.5,0.5)






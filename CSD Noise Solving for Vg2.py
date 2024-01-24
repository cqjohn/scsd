#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 14:19:39 2024

@author: genah
"""


#import libraries 
import matplotlib.pyplot as plt
import numpy as np 
import random

#define constants 
e = 1 #normalised everything so not 1.6e-19
pd = 50000 #is the number of 'points displayed' along 1 axis
pd_l = 6000 #same as pd but for transition lines as they may requrie less visible points 
pd_r = int(50000 / pd_l)  

#charging energies and coupling energy 
#Equations defined in Wiel et al 2002. 
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

#Solve Mu1 and Mu2 for Vg2 under condition that µ1 = µ2 
def Vg2(Vg1,N1,N2,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2):
    sumVg2 = (M1 - 0.5)*Ec1 + (M2 - N1)*Ecm -(N2-0.5)*Ec2 +(np.abs(1/e)*Cg1*Vg1*(Ecm-Ec1) + del_mu1 - del_mu2)
    ampVg2 = -np.abs(e)/ (Cg2*(Ec2-Ecm))
    
    vg2 = ampVg2 * sumVg2
    return vg2


#Solve Mu1 for Vg2 under condition that µ1 = 0 
def Vg2_mu1(Vg1,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1):
    sumVg2 = (M1 - 0.5)*Ec1 + (M2*Ecm) - (np.abs(1/e)*Cg1*Vg1*Ec1) + del_mu1
    ampVg2 = np.abs(e) / (Cg2*Ecm)
    return ampVg2 * sumVg2
    
#Solve Mu2 for Vg2 under condition that µ2 = 0 
def Vg2_mu2(Vg1,N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu2):
    sumVg2 = (N2 - 0.5)*Ec2 + N1*Ecm - (np.abs(1/e)*Cg1*Vg1*Ecm) + del_mu2 
    ampVg2 = np.abs(e) / (Cg2*Ec2)
    return sumVg2*ampVg2
    




#Solve for Mu1 given Vg1 and Vg2
def mu1_func(N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2,del_mu1):
    val = (N1 - 0.5)*Ec1 + N2*Ecm - (1/np.abs(e))*(Cg1*Vg1*Ec1 + Cg2*Vg2*Ecm) + del_mu1
    return val
#Solve for Mu2 given Vg1 and Vg2
def mu2_func(N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2,del_mu2):
    val = (N2 - 0.5)*Ec2 + N1*Ecm - (1/np.abs(e))*(Cg1*Vg1*Ecm + Cg2*Vg2*Ec2) + del_mu2
    return val


def plotChargeStab(CL,CR,Cm,Cg1,Cg2):
    
    #Eqs given in Wiel et al 2002.
    C1 = CL + Cm + Cg1
    C2 = CR + Cm + Cg2
    
    Ec1 = Ec1func(C1,C2,Cm)
    Ec2 = Ec2func(C1,C2,Cm)
    Ecm = Ecmfunc(C1,C2,Cm)
    
    #k is the mu1 electron occupancy 
    #j is the mu2 electron occupancy 
    for k in range(0,4):
        for j in range(0,4):
            
            #generate all Vg1 values
            Vg1 = np.linspace(0,10,pd)
            
            
            del_mu1 = [random.uniform(0, 1)*0.04 for _ in range(pd)]
            del_mu2 = [random.uniform(0, 1)*0.04 for _ in range(pd)]
            #del_mu1 = [0]*pd
            #del_mu2 = [0]*pd
            
            
            
            #1. FIND TRIPLE POINTS ------------------------------------
        
            #Electron type triple point  
  # Condition is Mu1(M1, M2) = Mu2(N1,N2) = where M1,M2,N1,N2 equals....
            M1e, M2e = k+1, j  #mu_1
            N1e, N2e = k, j+1  #mu_2
            
            #Hole type triple point 
  # Condition is Mu1(M1, M2) = Mu2(N1,N2) = where M1,M2,N1,N2 equals....
            M1h, M2h = k+1, j+1  #mu_1
            N1h, N2h = k+1, j+1  #mu_2
            
            #Find Vg2 for electron triple point 
            Vg2_e = Vg2(Vg1,k,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            #Find Vg2 for hole triple point  
            Vg2_h = Vg2(Vg1,N1h,N2h,M1h,M2h,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            
            #Find Mu1 and Mu2 for electron triple point
            #(Just to check if they are actually equal)
            mu1_e = mu1_func(M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2_e,del_mu1)
            mu2_e = mu2_func(N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2_e,del_mu2)
            
            #currently Vg2_e and Vg2_h hold all Mu_1 = Mu_2
            #however missing condition Mu_1 = Mu_2 = 0 hence must apply
            plot_Vg1_e= [] 
            plot_Vg2_e = []
            for i in range(1,len(mu1_e)):
                #print(mu1_e[i],mu2_e[i],Vg1[i],Vg2_e[i])
                if np.round(mu1_e[i],3) == 0:
                    #only plot the Vg1 and Vg2 that have a mu1/2=0 
                    plot_Vg1_e.append(Vg1[i])
                    plot_Vg2_e.append(Vg2_e[i])
                    break
                    
            #Find Mu1 and Mu2 for hole triple point
            #(Just to check if they are actually equal)        
            mu1_h = mu1_func(M1h,M2h,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2_h,del_mu1)
            mu2_h = mu2_func(N1h,N2h,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2_h,del_mu2)
             
            #same as above but for holes
            plot_Vg1_h = [] 
            plot_Vg2_h = []
           
            for i in range(1,len(mu1_h)):
                #print(mu1_h[i],mu2_h[i],Vg1[i],Vg2_h[i])
                if np.round(mu1_h[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_h.append(Vg1[i])
                     plot_Vg2_h.append(Vg2_h[i])
                     break
            
            #plot triple points         
            plt.plot(plot_Vg1_e,plot_Vg2_e,'.',c='black')
            plt.plot(plot_Vg1_h,plot_Vg2_h,'.',c='black')
            
            #2. FIND TRANSITION LINES ------------------------------------
            
            #2.1 CONDITION Mu_1 = MU_2 != 0 
            #Need to find the boundaries of the lines in terms of Vg1 
            lb = np.where(Vg1 == plot_Vg1_e[0])[0][0] #lower bound Vg1 
            ub = np.where(Vg1 == plot_Vg1_h[0])[0][0] #upper bound Vg1 
            #plot Vg2_e in this range, do not impose condition that mu1/2 = 0 
            plt.plot(Vg1[lb:ub:pd_r ],Vg2_e[lb:ub:pd_r],color='black')
            
           
            #2.2 CONDITION Mu_1 = 0 and Mu_2 = 0 
            #FInd Vg2 for condition
            Vg2_mu1_zero = Vg2_mu1(Vg1,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1)
            Vg2_mu2_zero = Vg2_mu2(Vg1,N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu2)
            
            #Find the adjacent hole triple points connected to the transition line ~~~~~~~~~~
            #Necessary to sefine boundaries
            M1hl, M2hl = k, j+1  #mu_1
            N1hl, N2hl = k, j+1  #mu_2
            
            M1hr, M2hr = k+1, j  #mu_1
            N1hr, N2hr = k+1, j  #mu_2
            
            #Additional lines 
            Vg2_hl = Vg2(Vg1,N1hl,N2hl,M1hl,M2hl,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            Vg2_hr = Vg2(Vg1,N1hr,N2hr,M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            
            plot_Vg1_hl = [] 
            plot_Vg2_hl = []
           
            mu1_hl = mu1_func(M1hl,M2hl,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2_hl,del_mu1)
            for i in range(1,len(mu1_hl)):
                if np.round(mu1_hl[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_hl.append(Vg1[i])
                     plot_Vg2_hl.append(Vg2_hl[i])
                     break
            
            plot_Vg1_hr = [] 
            plot_Vg2_hr = []        
            mu1_hr = mu1_func(M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2_hr,del_mu1)
            mu2_hr = mu2_func(M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2_hr,del_mu2)
            for i in range(1,len(mu1_hr)):
                if np.round(mu1_hr[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_hr.append(Vg1[i])
                     plot_Vg2_hr.append(Vg2_hr[i])
                     break
             #Adjacent hole triple points found so now we move on to define boundaries~~~ 
                     
            
          
            #Define upper and lower bounds of Vg1 for each transition line
            if not plot_Vg1_hl:
                print("Non found")
                lbl = 0 
                ubl = np.where(Vg1 == plot_Vg1_e[0])[0][0]
            else:
                lbl = np.where(Vg1 == plot_Vg1_hl[0])[0][0] 
                ubl = np.where(Vg1 == plot_Vg1_e[0])[0][0]
           
            if not plot_Vg1_hr: 
                print("Non found")
                lbr = np.where(Vg1 == plot_Vg1_e[0])[0][0]
                ubr = pd
            else: 
                lbr = np.where(Vg1 == plot_Vg1_e[0])[0][0]
                ubr = np.where(Vg1 == plot_Vg1_hr[0])[0][0]
            
        
            #plot adjacent hole triple points
            plt.plot(plot_Vg1_hl,plot_Vg2_hl,'.',color='white',  markeredgecolor='black')
            plt.plot(plot_Vg1_hr,plot_Vg2_hr,'.',color='white',  markeredgecolor='black')
            #plot transition lines with limits / boundaries 
            plt.plot(Vg1[lbr:ubr:pd_r ],Vg2_mu1_zero[lbr:ubr:pd_r ],c='black')
            plt.plot(Vg1[lbl:ubl:pd_r ],Vg2_mu2_zero[lbl:ubl:pd_r ],c='black')
           
            
           
            
           
            plt.xlim(0,8)
            plt.ylim(0,8)
            plt.title(f"Charge Stability Diagram")
            plt.xlabel("Vg1")
            plt.ylabel("Vg2") 
           
            
           
           
            
            
plotChargeStab(0.6,0.3,0.2,0.5,0.5)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 10:50:37 2024

@author: genah
"""


#import libraries 
import matplotlib.pyplot as plt
import numpy as np 
import random

#define constants 
e = 1 #normalised everything so not 1.6e-19
pd = 50000 #is the number of 'points displayed' along 1 axis
pd_l = 2000 #same as pd but for transition lines as they may requrie less visible points 
pd_r = int(50000 / pd_l)  



#Wiel et al. 2002 EQUATIONS ---------------------------------------------------



#SET 1: charging energies and coupling energy 
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



#SET 2: Solving for Vg2 given Vg1
#Solve µ1 and µ2 for Vg2 under condition that µ1 = µ2 
def Vg2(Vg1,N1,N2,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2):
    sumVg2 = (M1 - 0.5)*Ec1 + (M2 - N1)*Ecm -(N2-0.5)*Ec2 +(np.abs(1/e)*Cg1*Vg1*(Ecm-Ec1) + del_mu1 - del_mu2)
    ampVg2 = -np.abs(e)/ (Cg2*(Ec2-Ecm))
    vg2 = ampVg2 * sumVg2
    return vg2

#Solve µ1 for Vg2 under condition that µ1 = 0 
def Vg2_mu1(Vg1,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1):
    sumVg2 = (M1 - 0.5)*Ec1 + (M2*Ecm) - (np.abs(1/e)*Cg1*Vg1*Ec1) + del_mu1
    ampVg2 = np.abs(e) / (Cg2*Ecm)
    return ampVg2 * sumVg2
    
#Solve µ2 for Vg2 under condition that µ2 = 0 
def Vg2_mu2(Vg1,N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu2):
    sumVg2 = (N2 - 0.5)*Ec2 + N1*Ecm - (np.abs(1/e)*Cg1*Vg1*Ecm) + del_mu2 
    ampVg2 = np.abs(e) / (Cg2*Ec2)
    return sumVg2*ampVg2
    


#SET 3: Solving for Vg1 given Vg2
#Solve µ1 and µ2 for Vg1 under condition that µ1 = µ2 
def Vg1(Vg2,N1,N2,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2):
    sumVg1 = (M1 - 0.5)*Ec1 + (M2 - N1)*Ecm -(N2-0.5)*Ec2 +(np.abs(1/e)*Cg2*Vg2*(Ec2-Ecm) + del_mu1 - del_mu2)
    ampVg1 = -np.abs(e)/ (Cg1*(Ecm-Ec1))
    vg1 = ampVg1 * sumVg1
    return vg1

#Solve µ1 for Vg1 under condition that µ1 = 0 
def Vg1_mu1(Vg2,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1):
    sumVg1 = (M1 - 0.5)*Ec1 + (M2*Ecm) - (np.abs(1/e)*Cg2*Vg2*Ecm) + del_mu1
    ampVg1 = np.abs(e) / (Cg1*Ec1)
    return ampVg1 * sumVg1
  
#Solve µ2 for Vg1 under condition that µ2 = 0 
def Vg1_mu2(Vg2,N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu2):
    sumVg1 = (N2 - 0.5)*Ec2 + N1*Ecm - (np.abs(1/e)*Cg2*Vg2*Ec2) + del_mu2 
    ampVg1 = np.abs(e) / (Cg1*Ecm)
    return sumVg1*ampVg1
    


#SET 4: µ1 and µ2 given Vg1 and Vg2
#Solve for Mu1 given Vg1 and Vg2
def mu1_func(N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2,del_mu1):
    val = (N1 - 0.5)*Ec1 + N2*Ecm - (1/np.abs(e))*(Cg1*Vg1*Ec1 + Cg2*Vg2*Ecm) + del_mu1
    return val
#Solve for Mu2 given Vg1 and Vg2
def mu2_func(N1,N2,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1,Vg2,del_mu2):
    val = (N2 - 0.5)*Ec2 + N1*Ecm - (1/np.abs(e))*(Cg1*Vg1*Ecm + Cg2*Vg2*Ec2) + del_mu2
    return val



#SET 5: Solving for V2 and V1 with known µ1 and µ2
def Vg2_fmu(N1,N2,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,mu1,mu2,del_mu1, del_mu2):
    amp = Ecm/Ec1 
    sum1 = ((M1-0.5)*Ec1) + (M2*Ecm) + del_mu1 - mu1
    
    ampVg2 = (Ec1 / ((Ec1*Ec2) - (Ecm**2))) * (np.abs(e)/Cg2)
    sumVg2 = (N2-0.5)*Ec2 + (N1*Ecm) + del_mu2 - mu2 - amp*sum1
    return sumVg2 * ampVg2
  
     
    
def Vg1_fmu(N1,N2,M1,M2,Ec1,Ec2,Ecm,Cg1,Cg2,mu1,mu2,del_mu1, del_mu2,Vg2):  
    ampVg1 = np.abs(e) / (Cg1*Ec1)
    sumVg1 = ((M1 - 0.5)*Ec1) + (M2*Ecm) - ((1/np.abs(e))*Cg2*Vg2*Ecm) + del_mu1 - mu1
    return ampVg1 * sumVg1
#------------------------------------------------------------------------------


#Ploting the charge stability diagram (CSD)
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
            
            #generate all Vg1/Vg2 values
            Vg2_set = np.linspace(0,10,pd)
            Vg1_set = np.linspace(0,10,pd)
            
            #calculate no noise CSD
            del_mu1 = [0]*pd
            del_mu2 = [0]*pd
        
            
            
#NOISELESS CSD: used to find the values that we want to work with

#(SET VG2) STEP 1: FIND TRIPLE POINTS SET Vg2
        
            #Electron type triple point  
            # Condition is Mu1(M1, M2) = Mu2(N1,N2) = where M1,M2,N1,N2 equals....
            M1e, M2e = k+1, j  #mu_1
            N1e, N2e = k, j+1  #mu_2
            
            #Hole type triple point 
            # Condition is Mu1(M1, M2) = Mu2(N1,N2) = where M1,M2,N1,N2 equals....
            M1h, M2h = k+1, j+1  #mu_1
            N1h, N2h = k+1, j+1  #mu_2
            
            #Find Vg1 for electron triple point 
            Vg1_e = Vg1(Vg2_set,N1e,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
          
            #Find Vg2 for hole triple point  
            Vg1_h = Vg1(Vg2_set,N1h,N2h,M1h,M2h,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            
            #Find Mu1 and Mu2 for electron triple point
            #(Just to check if they are actually equal)
            mu1_e_svg2 = mu1_func(M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_e,Vg2_set,del_mu1)
            mu2_e_svg2 = mu2_func(N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_e,Vg2_set,del_mu2)
           
            #currently Vg2_e and Vg2_h hold all Mu_1 = Mu_2
            #however missing condition Mu_1 = Mu_2 = 0 hence must apply
            plot_Vg1_e_svg2= [] 
            plot_Vg2_e_svg2 = []
            for i in range(1,len(mu1_e_svg2)):
                #print(mu1_e[i],mu2_e[i],Vg1[i],Vg2_e[i])
                if np.round(mu1_e_svg2[i],3) == 0:
                    #only plot the Vg1 and Vg2 that have a mu1/2=0 
                    plot_Vg1_e_svg2.append(Vg1_e[i])
                    plot_Vg2_e_svg2.append(Vg2_set[i])
                    break
                    
            #Find Mu1 and Mu2 for hole triple point
            #(Just to check if they are actually equal)        
            mu1_h_svg2 = mu1_func(M1h,M2h,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_h,Vg2_set,del_mu1)
            mu2_h_svg2 = mu2_func(N1h,N2h,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_h,Vg2_set,del_mu2)
             
            #same as above but for holes
            plot_Vg1_h_svg2 = [] 
            plot_Vg2_h_svg2 = []
           
            for i in range(1,len(mu1_h_svg2)):
                #print(mu1_h[i],mu2_h[i],Vg1[i],Vg2_h[i])
                if np.round(mu1_h_svg2[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_h_svg2.append(Vg1_h[i])
                     plot_Vg2_h_svg2.append(Vg2_set[i])
                     break
            

#(SET VG1) STEP 1: FIND TRIPLE POINTS SET Vg1
            #Electron type triple point  
            # Condition is Mu1(M1, M2) = Mu2(N1,N2) = where M1,M2,N1,N2 equals....
            M1e, M2e = k+1, j  #mu_1
            N1e, N2e = k, j+1  #mu_2
            
            #Hole type triple point 
            # Condition is Mu1(M1, M2) = Mu2(N1,N2) = where M1,M2,N1,N2 equals....
            M1h, M2h = k+1, j+1  #mu_1
            N1h, N2h = k+1, j+1  #mu_2
            
            #Find Vg2 for electron triple point 
            
            Vg2_e = Vg2(Vg1_set,k,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            #Find Vg2 for hole triple point  
            Vg2_h = Vg2(Vg1_set,N1h,N2h,M1h,M2h,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            
            #Find Mu1 and Mu2 for electron triple point
            #(Just to check if they are actually equal)
            mu1_e_svg1 = mu1_func(M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_set,Vg2_e,del_mu1)
            mu2_e_svg1 = mu2_func(N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_set,Vg2_e,del_mu2)
            
            #currently Vg2_e and Vg2_h hold all Mu_1 = Mu_2
            #however missing condition Mu_1 = Mu_2 = 0 hence must apply
            plot_Vg1_e_svg1 = [] 
            plot_Vg2_e_svg1 = []
            for i in range(1,len(mu1_e_svg1)):
                if np.round(mu1_e_svg1[i],3) == 0:
                    #only plot the Vg1 and Vg2 that have a mu1/2=0 
                    plot_Vg1_e_svg1.append(Vg1_set[i])
                    plot_Vg2_e_svg1.append(Vg2_e[i])
                    break
                    
            #Find Mu1 and Mu2 for hole triple point
            #(Just to check if they are actually equal)        
            mu1_h_svg1 = mu1_func(M1h,M2h,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_set,Vg2_h,del_mu1)
            mu2_h_svg1 = mu2_func(N1h,N2h,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_set,Vg2_h,del_mu2)
             
            #same as above but for holes
            plot_Vg1_h_svg1 = [] 
            plot_Vg2_h_svg1 = []
           
            for i in range(1,len(mu1_h_svg1)):
                if np.round(mu1_h_svg1[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_h_svg1.append(Vg1_set[i])
                     plot_Vg2_h_svg1.append(Vg2_h[i])
                     break
           

#(SET VG2 STEP 2.1)   
            #2.1 CONDITION Mu_1 = MU_2 != 0  SET Vg2
            #Need to find the boundaries of the lines in terms of Vg1 
            lb_svg2 = np.where(Vg2_set == plot_Vg2_e_svg2[0])[0][0] #lower bound Vg1 
            ub_svg2 = np.where(Vg2_set == plot_Vg2_h_svg2[0])[0][0] #upper bound Vg1 
    

#(SET VG1 STEP 2.1)  
            #2.1 CONDITION Mu_1 = MU_2 != 0 
            #Need to find the boundaries of the lines in terms of Vg1 
            lb_svg1 = np.where(Vg1_set == plot_Vg1_e_svg1[0])[0][0] #lower bound Vg1 
            ub_svg1 = np.where(Vg1_set == plot_Vg1_h_svg1[0])[0][0] #upper bound Vg1 
           
            

#(SET VG2 STEP 2.2)  
            #2.2 CONDITION Mu_1 = 0 and Mu_2 = 0 set Vg2 
            
            #FInd Vg2 for condition
            Vg1_mu1_zero = Vg1_mu1(Vg2_set,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1)
            Vg1_mu2_zero = Vg1_mu2(Vg2_set,N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu2)
            
            #Find the adjacent hole triple points connected to the transition line 
            #Necessary to sefine boundaries
            M1hl, M2hl = k, j+1  #mu_1
            N1hl, N2hl = k, j+1  #mu_2
            
            M1hr, M2hr = k+1, j  #mu_1
            N1hr, N2hr = k+1, j  #mu_2
            
            #Additional lines 
            Vg1_hl = Vg1(Vg2_set,N1hl,N2hl,M1hl,M2hl,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            Vg1_hr = Vg1(Vg2_set,N1hr,N2hr,M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            
            plot_Vg1_hl_svg2 = [] 
            plot_Vg2_hl_svg2 = []
           
            mu1_hl_svg2 = mu1_func(M1hl,M2hl,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_hl,Vg2_set,del_mu1)
            for i in range(1,len(mu1_hl_svg2)):
                if np.round(mu1_hl_svg2[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_hl_svg2.append(Vg1_hl[i])
                     plot_Vg2_hl_svg2.append(Vg2_set[i])
                     break
            
            plot_Vg1_hr_svg2 = [] 
            plot_Vg2_hr_svg2 = []        
            mu1_hr_svg2 = mu1_func(M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_hr,Vg2_set,del_mu1)
            mu2_hr_svg2 = mu2_func(M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_hr,Vg2_set,del_mu2)
            for i in range(1,len(mu1_hr_svg2)):
                if np.round(mu1_hr_svg2[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_hr_svg2.append(Vg1_hr[i])
                     plot_Vg2_hr_svg2.append(Vg2_set[i])
                     break
            #Adjacent hole triple points found so now we move on to define boundaries 
            

        
            #Define upper and lower bounds of Vg1 for each transition line
            #left 
            if not plot_Vg1_hl_svg2:
                lbl_svg2 = np.where(Vg2_set == plot_Vg2_e_svg2[0])[0][0] 
                ubl_svg2 = pd
            else:
                lbl_svg2 = np.where(Vg2_set == plot_Vg2_e_svg2[0])[0][0] 
                ubl_svg2 = np.where(Vg2_set == plot_Vg2_hl_svg2[0])[0][0]
           
            if not plot_Vg1_hr_svg2: 
                lbr_svg2 = 0
                ubr_svg2 = np.where(Vg2_set == plot_Vg2_e_svg2[0])[0][0]
            else: 
                lbr_svg2 = np.where(Vg2_set == plot_Vg2_hr_svg2[0])[0][0]
                ubr_svg2 = np.where(Vg2_set == plot_Vg2_e_svg2[0])[0][0]
            

            
            

#(SET VG1 STEP 2.2)  
            
            #2.2 CONDITION Mu_1 = 0 and Mu_2 = 0 set Vg1 
            
            
            #FInd Vg2 for condition
            Vg2_mu1_zero = Vg2_mu1(Vg1_set,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1)
            Vg2_mu2_zero = Vg2_mu2(Vg1_set,N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu2)
            
            #Find the adjacent hole triple points connected to the transition line 
            #Necessary to sefine boundaries
            M1hl, M2hl = k, j+1  #mu_1
            N1hl, N2hl = k, j+1  #mu_2
            
            M1hr, M2hr = k+1, j  #mu_1
            N1hr, N2hr = k+1, j  #mu_2
            
            #Additional lines 
            Vg2_hl = Vg2(Vg1_set,N1hl,N2hl,M1hl,M2hl,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            Vg2_hr = Vg2(Vg1_set,N1hr,N2hr,M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,del_mu1,del_mu2)
            
            plot_Vg1_hl_svg1 = [] 
            plot_Vg2_hl_svg1 = []
           
            mu1_hl_svg1 = mu1_func(M1hl,M2hl,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_set,Vg2_hl,del_mu1)
            for i in range(1,len(mu1_hl_svg1)):
                if np.round(mu1_hl_svg1[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_hl_svg1.append(Vg1_set[i])
                     plot_Vg2_hl_svg1.append(Vg2_hl[i])
                     break
            
            plot_Vg1_hr_svg1 = [] 
            plot_Vg2_hr_svg1 = []        
            mu1_hr_svg1 = mu1_func(M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_set,Vg2_hr,del_mu1)
            mu2_hr_svg1 = mu2_func(M1hr,M2hr,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_set,Vg2_hr,del_mu2)
            for i in range(1,len(mu1_hr_svg1)):
                if np.round(mu1_hr_svg1[i],3) == 0:
                     #print(mu1_h[i])
                     plot_Vg1_hr_svg1.append(Vg1_set[i])
                     plot_Vg2_hr_svg1.append(Vg2_hr[i])
                     break
             #Adjacent hole triple points found so now we move on to define boundaries
                     
          
            #Define upper and lower bounds of Vg1 for each transition line
            if not plot_Vg1_hl_svg1:
                lbl_svg1 = 0 
                ubl_svg1 = np.where(Vg1_set == plot_Vg1_e_svg1[0])[0][0]
            else:
                lbl_svg1 = np.where(Vg1_set == plot_Vg1_hl_svg1[0])[0][0] 
                ubl_svg1 = np.where(Vg1_set == plot_Vg1_e_svg1[0])[0][0]
           
            if not plot_Vg1_hr_svg1: 
                lbr_svg1 = np.where(Vg1_set == plot_Vg1_e_svg1[0])[0][0]
                ubr_svg1 = pd
            else: 
                lbr_svg1 = np.where(Vg1_set == plot_Vg1_e_svg1[0])[0][0]
                ubr_svg1 = np.where(Vg1_set == plot_Vg1_hr_svg1[0])[0][0]
             
#Plot noiseless data 
            '''
            #con1 is µ1 = µ2 
            plt.plot(Vg1_e[lb_svg2:ub_svg2:pd_r ],Vg2_set[lb_svg2:ub_svg2:pd_r],'.',color='black')
            #con2 is µ1 = 0 
            plt.plot(Vg1_mu1_zero[lbr_svg2:ubr_svg2:pd_r ],Vg2_set[lbr_svg2:ubr_svg2:pd_r ],'.',c='black')
            #con3 ia µ2 = 0 
            plt.plot(Vg1_set[lbl_svg1:ubl_svg1:pd_r ],Vg2_mu2_zero[lbl_svg1:ubl_svg1:pd_r ],'.',c='black')
            '''
            
            
#NOISY DATA 
            
#STEP 1     #Get noiseless µ's 
            #con1 is µ1 = µ2 
            len_con1 = len(Vg2_set[lb_svg2:ub_svg2:pd_r])
            del_mu1 = np.array([0]*len_con1)
            del_mu2 = np.array([0]*len_con1)
            mu1_con1 = mu1_func(M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_e[lb_svg2:ub_svg2:pd_r ],Vg2_set[lb_svg2:ub_svg2:pd_r],del_mu1 )
            mu2_con1 = mu2_func(N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_e[lb_svg2:ub_svg2:pd_r ],Vg2_set[lb_svg2:ub_svg2:pd_r],del_mu2 )
            
            #con2 is µ1 = 0 
            len_con2 = len(Vg2_set[lbr_svg2:ubr_svg2:pd_r ])
            del_mu1 = np.array([0]*len_con2)
            del_mu2 = np.array([0]*len_con2)
            mu1_con2 = mu1_func(M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_mu1_zero[lbr_svg2:ubr_svg2:pd_r ],Vg2_set[lbr_svg2:ubr_svg2:pd_r ], del_mu1)
            mu2_con2 = mu2_func(N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2,Vg1_mu1_zero[lbr_svg2:ubr_svg2:pd_r ],Vg2_set[lbr_svg2:ubr_svg2:pd_r ], del_mu2)
            
            #con3 ia µ2 = 0 
            len_con3 = len(Vg2_mu2_zero[lbl_svg1:ubl_svg1:pd_r ])  
            del_mu1 = np.array([0]*len_con3)
            del_mu2 = np.array([0]*len_con3)                  
            mu1_con3 = mu1_func(M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2, Vg1_set[lbl_svg1:ubl_svg1:pd_r ],Vg2_mu2_zero[lbl_svg1:ubl_svg1:pd_r ], del_mu1 )
            mu2_con3 = mu2_func(N1e,N2e,Ec1,Ec2,Ecm,Cg1,Cg2, Vg1_set[lbl_svg1:ubl_svg1:pd_r ],Vg2_mu2_zero[lbl_svg1:ubl_svg1:pd_r ],del_mu2)
            
            
            
#STEP 2:    #Find Vg1 and Vg2 with known µ1, µ2, ∆µ1, ∆µ2 for all three conditions 
            del_mu1 = np.array([random.uniform(-1, 1)*0.04 for _ in range(pd)])
            del_mu2 = np.array([random.uniform(-1, 1)*0.04 for _ in range(pd)])
            
            #con1 is µ1 = µ2
            Vg2_con1 = Vg2_fmu(N1e,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,mu1_con1,mu2_con1, del_mu1[lb_svg2:ub_svg2:pd_r], del_mu2[lb_svg2:ub_svg2:pd_r])
            Vg1_con1 = Vg1_fmu(N1e,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,mu1_con1,mu2_con1, del_mu1[lb_svg2:ub_svg2:pd_r], del_mu2[lb_svg2:ub_svg2:pd_r],Vg2_con1)
            
            #con2 is µ1 = 0 
            Vg2_con2 = Vg2_fmu(N1e,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,mu1_con2,mu2_con2, del_mu1[lbr_svg2:ubr_svg2:pd_r ], del_mu2[lbr_svg2:ubr_svg2:pd_r ])
            Vg1_con2 = Vg1_fmu(N1e,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,mu1_con2,mu2_con2, del_mu1[lbr_svg2:ubr_svg2:pd_r ], del_mu2[lbr_svg2:ubr_svg2:pd_r ],Vg2_con2)
            
            #con3 ia µ2 = 0 
            Vg2_con3 = Vg2_fmu(N1e,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,mu1_con3,mu2_con3, del_mu1[lbl_svg1:ubl_svg1:pd_r ], del_mu2[lbl_svg1:ubl_svg1:pd_r ])
            Vg1_con3 = Vg1_fmu(N1e,N2e,M1e,M2e,Ec1,Ec2,Ecm,Cg1,Cg2,mu1_con3,mu2_con3, del_mu1[lbl_svg1:ubl_svg1:pd_r ], del_mu2[lbl_svg1:ubl_svg1:pd_r ],Vg2_con3)
            
            
#PLOT NOISY DATA 
            plt.plot(Vg1_con1,Vg2_con1,'.',color= 'purple')
            plt.plot(Vg1_con2,Vg2_con2,'.',color='purple')
            plt.plot(Vg1_con3,Vg2_con3,'.',color= 'purple')
            
            
            plt.xlim(0,8)
            plt.ylim(0,8)
            plt.title(f"Charge Stability Diagram")
            plt.xlabel("Vg1")
            plt.ylabel("Vg2") 
            
            
            
            
            
            
            
            
            
plotChargeStab(0.6,0.3,0.0000001,0.5,0.5)

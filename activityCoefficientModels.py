#------------------------------------------------------------------------------------------------------------
#
#
#                 __  __  ____  _____  ______ _       _____   ______ ____  _____           
#                |  \/  |/ __ \|  __ \|  ____| |     / ____| |  ____/ __ \|  __ \          
#                | \  / | |  | | |  | | |__  | |    | (___   | |__ | |  | | |__) |         
#                | |\/| | |  | | |  | |  __| | |     \___ \  |  __|| |  | |  _  /          
#                | |  | | |__| | |__| | |____| |____ ____) | | |   | |__| | | \ \          
#                |_|  |_|\____/|_____/|______|______|_____/  |_|    \____/|_|  \_\         
#                           _____ _______ _______      _______ _________     __             
#                     /\   / ____|__   __|_   _\ \    / /_   _|__   __\ \   / /             
#                    /  \ | |       | |    | |  \ \  / /  | |    | |   \ \_/ /              
#                   / /\ \| |       | |    | |   \ \/ /   | |    | |    \   /               
#                  / ____ \ |____   | |   _| |_   \  /   _| |_   | |     | |                
#                 /_/    \_\_____|  |_|  |_____|   \/   |_____|  |_|     |_|                
#             _____ ____  ______ ______ ______ _____ _____ _____ ______ _   _ _______ 
#            / ____/ __ \|  ____|  ____|  ____|_   _/ ____|_   _|  ____| \ | |__   __|
#           | |   | |  | | |__  | |__  | |__    | || |      | | | |__  |  \| |  | |   
#           | |   | |  | |  __| |  __| |  __|   | || |      | | |  __| | . ` |  | |   
#           | |___| |__| | |____| |    | |     _| || |____ _| |_| |____| |\  |  | |   
#            \_____\____/|______|_|    |_|    |_____\_____|_____|______|_| \_|  |_|   
#                                                                                                                                                      
#                                                                                                                                            
#                                                                                                                                            
#------------------------------------------------------------------------------------------------------------
#
#   Technische Universiteit Delft - TUDelft (2022)
#
#   Master of Science in Chemical Engineering
#
#   This code was developed and tested by ELIA FERRETTI
#
#   You can redistribute the code and/or modify it
#   Whenever the code is used to produce any publication or document,
#   reference to this work and author should be reported
#   No warranty of fitness for a particular purpose is offered
#   The user must assume the entire risk of using this code
#
#------------------------------------------------------------------------------------------------------------
#
#   student name:       ELIA FERRETTI
#   e-mail:             eliaferretti@outlook.it
#   LinkedIn:           https://www.linkedin.com/in/elia-ferretti/
#   GitHub:             https://github.com/eliaferretti
#   Instagram:          @eliaferretti
#
#------------------------------------------------------------------------------------------------------------


#import needed library
import numpy as np

def activityCoefficient_WILSON(i,x,W):
    N = len(x)
    alpha = 0
    gamma = 0
    for j in range(0,N):
        alpha = alpha + x[j]*W[i,j]
    
    for k in range(0,N):
        beta = 0
        for j in range(0,N):
            beta = beta + x[j]*W[k,j]
        gamma = gamma + x[k]*W[k,i]/beta
    
    gamma = np.exp( - np.log(alpha) + 1 - gamma )

    return gamma

def activityCoefficient_NRTL(i,x,G,tau):
    
    N = len(x)
    alpha = 0
    beta = 0
    phi = 0
    
    for j in range(0,N):
        alpha = alpha + tau[j,i]*G[j,i]*x[j]
        beta = beta + G[j,i]*x[j]
    
    for j in range(0,N):
        gamma = 0; delta = 0; epsilon = 0
        for k in range(0,N):
            gamma = gamma + G[k,j]*x[k]
            delta= delta + x[k]*tau[k,j]*G[k,j]
            epsilon = epsilon + x[k]*G[k,j]
        
        phi = phi + x[j]*G[i,j]/gamma*( tau[i,j] - delta/epsilon )
            
    gamma = np.exp( alpha/beta + phi )
    
    return gamma
    
def activityCoefficient_UNIQUAC(i,x,tau,z,r,q):
    
    N = len(x)
    
    I = z/2*(r-q) - (r - 1)
    
    alpha = 0; beta = 0; gamma = 0; delta = 0; rho = 0
    for j in range(0,N):
        alpha += q[j]*x[j]
        beta += r[j]*x[j]
        gamma += x[j]*I[j]
          
    theta = q*x/alpha
    phi = r*x/beta
    
    for j in range(0,N):       
        epsilon = 0
        for k in range(0,N):
            epsilon += theta[k]*tau[k,j]
            
        delta += theta[j]*tau[j,i]
        rho += theta[j]*tau[i,j]/epsilon  
    
    logGamma_C = np.log(r[i]/beta) + z/2*q[i]*np.log(q[i]*beta/r[i]/alpha) + I[i] - r[i]/beta*gamma   
    logGamma_R = q[i]*( 1 - np.log(delta) - rho ) 
    
    return np.exp(logGamma_C + logGamma_R)

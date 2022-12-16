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

import numpy as np
import scipy
import matplotlib.pyplot as plt

spec = ["tert-butanol","toluene"]
method = ["Wilson", "NRTL", "UNIQUAC"]

#experimental data @ T = 35°C
A = np.array([[1028.5272,   491.8786,   0],
              [-3225.6147,  4711.5914,  0.0159],
              [25.9015,     316.6514,   0]])
T = 35 + 273.15

#data
vL = np.array([94.88, 106.85])
R = np.array([3.4528, 3.9228])
Q = np.array([3.128, 2.968])

#wilson parameter
W = np.ones((len(spec),len(spec)))
W[0,1] = vL[1]/vL[0]*np.exp(-A[0,0]/1.98721/T)
W[1,0] = vL[0]/vL[1]*np.exp(-A[0,1]/1.98721/T)

#NRTL parameter
tau = np.zeros((len(spec),len(spec)))
G = np.ones((len(spec),len(spec)))

tau[0,1] = A[1,0]/1.98721/T
tau[1,0] = A[1,1]/1.98721/T
G[0,1] = np.exp(-A[1,2]*tau[0,1])
G[1,0] = np.exp(-A[1,2]*tau[1,0])

#UNIQUAC parameter
z = 10
tauN = np.ones((len(spec),len(spec)))
tauN[0,1] = np.exp(-A[2,0]/1.98721/T)
tauN[1,0] = np.exp(-A[2,1]/1.98721/T)

#Anotine parameter
A = np.array([7.23159, 6.95087])
B = np.array([1107.06, 1342.31])
C = np.array([172.101, 219.187])


#------------------------------------------------------------------------------
#                           ---TESTING---
#------------------------------------------------------------------------------


n = 50
x = np.linspace(0,1,n)
X = np.transpose(np.array([x,1-x]))

p0 = np.power(10, A - B/(T+C-273.15))

y = np.zeros(n); y1 = np.zeros(n); y2 = np.zeros(n)
pp = np.zeros((n,3))

for i in range(0,n):
    eqn = lambda p: 1 - x[i]*p0[0]*activityCoefficient_WILSON(0,X[i],W)/p - (1-x[i])*p0[1]*activityCoefficient_WILSON(1,X[i],W)/p
    eqn2 = lambda p: 1 - x[i]*p0[0]*activityCoefficient_NRTL(0,X[i],G,tau)/p - (1-x[i])*p0[1]*activityCoefficient_NRTL(1,X[i],G,tau)/p
    eqn3 = lambda p: 1 - x[i]*p0[0]*activityCoefficient_UNIQUAC(0,X[i],tauN,z,R,Q)/p - (1-x[i])*p0[1]*activityCoefficient_UNIQUAC(1,X[i],tauN,z,R,Q)/p
    
    p = scipy.optimize.fsolve(eqn,p0[0])
    p1 = scipy.optimize.fsolve(eqn2,p0[0])
    p2 = scipy.optimize.fsolve(eqn3,p0[0])
    
    pp[i] = p[0], p1[0], p2[0]
    
    y[i] = x[i]*p0[0]*activityCoefficient_WILSON(0,X[i],W)/p
    y1[i] = x[i]*p0[0]*activityCoefficient_NRTL(0,X[i],G,tau)/p1
    y2[i] = x[i]*p0[0]*activityCoefficient_UNIQUAC(0,X[i],tauN,z,R,Q)/p2
    
plt.figure(1,figsize=(5,5))
plt.plot(x,y,"-b",x,y1,"-g",x,y2,"-r",x,x,":k")
plt.title("Equilibrium curve at T = 35 °C")
plt.xlabel("x [-], tert-butanol"); plt.ylabel("y [-], tert-butanol")
plt.legend(["WILSON","NRTL","UNIQUAC"])
plt.grid(); plt.show()
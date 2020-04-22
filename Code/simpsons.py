#!/usr/bin/env python3

"""

@file simpsons.py
@brief Simpson's Rule for integration of a function

This file contains a basic 1D setup for the Simpson's Rule for integration of a function. This
method relies on approximation of a given function in terms of Lagrange polynomials and
then calculating the area under the approximated function. 

==============================================================================

@author Sankarasubramanian Ragunathan
@date 14.04.2020

"""

# Importing packages required for numerical calculations and plotting
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')

# Function definition to perform Simpson's Integration
def simpsons(a,b,N) :
    
    # Grid spacing for the Uniform grid of points
    h = (b-a)/(N-1)            
    # Area obtained from integrating f(x) using Simpson's Rule                         
    area = 0                                            
    
    # Checking if the number of slices are even or odd (necessary for the implementation of Simpson's Rule)
    if (N%2 == 1) :
        # Looping through all the points on the grid
        for k in range(1,N-1,2) :
        	# Point x_k to evaluate the function value f(x_k)
            xk = a+k*h          
            # Calculating the area under the sin(x) curve using Simpson's Rule                        
            area += (h/3)*(np.sin(xk-h) + \
                           4*np.sin(xk) +  \
                           np.sin(xk+h))                
    else :
        # Looping through all the points on the grid
        for k in range(1,N-2,2) :
        	# Point x_k to evaluate the function value f(x_k)
            xk = a+k*h      
            # Calculating the area under the sin(x) curve using Simpson's Rule                            
            area += (h/3)*(np.sin(xk-h) + \
                           4*np.sin(xk) + \
                           np.sin(xk+h))     
        # Special treatment of the boundary term for even number of points in the grid           
        xk = a + (N-2)*h
        # Calculating the area under the sin(x) curve using Simpson's Rule                                
        area += (h/12)*(-np.sin(xk-h) + \
                        8*np.sin(xk) + \
                        5*np.sin(xk+h))                 
    return area

# Function defined to calculate the Analytical error and Theoretical error
def errorFunc(a,b,N,area,errOpt) :
	# Grid spacing for the uniform grid
    h = (b-a)/(N-1)
    # Actual solution obtained from integrating the function                                     
    actualSoln = np.cos(a) - np.cos(b)                  
    if (errOpt == 0) :
    	# Analytical error calculation
        error = np.fabs(actualSoln-area)                
    elif (errOpt == 1) :
    	# Theoretical discretization error
        constFac = 10
        error = (1/180)*(b-a)*(h**4)*constFac           
    return (h,error)

""" Main Program """

# No. of terms to use in the generation of array of log-spaced grid points
nPoints = 20       
# Array containing the integral calculated using Simpson's Rule for different no. of grid points                                                  
area = np.zeros(nPoints)
# log-spaced values for the no. of gridpoints in the uniform mesh                                             
N = np.logspace(1,6,num=nPoints,endpoint=True,dtype=np.int64)        
for index in range(0,nPoints) :
    area[index] = simpsons(0, np.pi/2, N[index])

(h,analyticalerror) = errorFunc(0,np.pi/2,N,area,0)    
(h,theoreticalerror) = errorFunc(0,np.pi/2,N,area,1)

# Plotting the error vs. grid spacing (log-log plot)

plt.loglog(h,analyticalerror,label=r"$|S-S_{S}(h)|$",linewidth=2,color="r",marker="o",markersize=8)
plt.loglog(h,theoreticalerror,label=r"$\mathcal{O}(h^4)$",linewidth=2,color="k",marker="s",markersize=6)
plt.grid(True,which="both",axis="both",linestyle=":")
plt.title(r"Simpson's Rule",fontsize=16)
plt.xlabel(r'$h = \frac{b-a}{N-1}$',fontsize=14)
plt.ylabel(r'Error',fontsize=14)
plt.legend(loc=2,fontsize=14)
plt.savefig('./simpsons.eps',format='eps')

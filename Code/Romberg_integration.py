"""

@file Romberg_integration.py
@brief Romberg integrator incorporating Neville scheme

This file contains a basic 1D setup for the Romberg integration scheme. This
method relies on a Richardson extrapolation of integrals calculated with the
trapezoidal rule on a series of successively refined grids.
The practical implementation relies on the Neville scheme.  Three test
functions are provided, of which the error is visualized over the number of 
levels in the Romberg integration procedure.

==============================================================================

@author Maximilian Kruse
@date 14.04.2020

"""

#============================= Import Libraries ==============================
import numpy as np
import matplotlib.pyplot as plt

 
#========================== Integration Functions ============================

#-----------------------------------------------------------------------------
"""
Calculate Trapezoidal integral

@param[in]: function handle, number of grid points, integration domain bounds
@param[out]: Value of the definite Trapezoidal integral
"""
def trapezoidal_integral(funcHandle, numPoints, lb, ub):
    # Calculate function values on grid
    discrSpacing = (ub-lb)/(numPoints-1)
    discrPoints = [lb+(ub-lb)*i/(numPoints-1) for i in range(numPoints)]
    funcVals = funcHandle(discrPoints)
    
    # Accumulate according to trapezoidal rule
    integralValue = 0;
    for i  in range(numPoints-1):
        integralValue += funcVals[i] + funcVals[i+1]
    integralValue *= discrSpacing/2;
    
    return integralValue

#-----------------------------------------------------------------------------
"""
Calculate Romberg integral

@param[in]: function handle, number levels, integration domain bounds
@param[out]: Value of the definite Romberg integral
Note: This function calls the trapezoidal_integral function for intrinsic
      calculations
"""
def romberg_integral(funcHandle, levels, lb, ub):
    # Define initial column of neville scheme from trapezoidal rule
    schemeArray = [trapezoidal_integral(funcHandle, pow(2,i)+1, lb, ub)
                   for i in range(levels+1)]
    
    # Successively compute new columns of the Neville scheme. This routine 
    # overwrites "non-diagonal" elements for reduced storage requirements
    for level in range(1, levels+1):
        for i in range (levels,level-1,-1):
            schemeArray[i] = schemeArray[i] + 1/(pow(4,level)-1)*\
                            (schemeArray[i]-schemeArray[i-1])
                             
    return np.asarray(schemeArray)


#============================== Test Functions ===============================

# Function handles
expFunc = lambda x:np.exp(x)
sineFunc = lambda x:np.power(np.sin(8*x),4)
sqrtFunc = lambda x:np.sqrt(x)

# Analytical solutions  
solExp = np.exp(1)-1
solSine = 3./4.*np.pi;
solSqrt = 2./3.


#============================= Main Computation ==============================

# Set levels for the test integrals
levelExp = 4
levelSine = 8
levelSqrt = 20

# Compute Romberg integrals for the test functions
xValsExp = np.arange(0, levelExp+1)
xValsSine = np.arange(0, levelSine+1)
xValsSqrt = np.arange(0, levelSqrt+1)
yValsExp = np.absolute(romberg_integral(expFunc, levelExp, 0, 1)-solExp)
yValsSine = np.absolute(romberg_integral(sineFunc, levelSine, 0, 2*np.pi)-solSine)
yValsSqrt = np.absolute(romberg_integral(sqrtFunc, levelSqrt, 0, 1)-solSqrt)

# Plot results
plt.title('Error of the Romberg integration scheme')
plt.xlabel(r'Level $k$')
plt.ylabel(r'Absolute Error $|\epsilon|$')
plt.xticks(np.arange(0, 21, step=5))
plt.grid('on')
plt.semilogy(xValsExp,yValsExp,label=r'$e^x$',marker='o',color='tab:blue',)
plt.semilogy(xValsSine,yValsSine,label=r'$sin^4(8x)$',marker='o',color='tab:orange')
plt.semilogy(xValsSqrt,yValsSqrt,label=r'$\sqrt{x}$',marker='o',color='tab:green')
plt.legend()
plt.show()


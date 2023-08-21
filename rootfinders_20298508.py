"""
MATH2050 CW1 rootfinders module

@author: Saanvi Mudholkar
"""


import numpy as np
import matplotlib.pyplot as plt

#Question 2
def bisection(f,a,b,Nmax):
    
    """
    Bisection Method: Returns a numpy array of the 
    sequence of approximations obtained by the bisection method.
    
    Parameters
    ----------
    f : function
        Input function for which the zero is to be found.
    a : real number
        Left side of interval.
    b : real number
        Right side of interval.
    Nmax : integer
        Number of iterations to be performed.
        
    Returns
    -------
    p_array : numpy.ndarray, shape (Nmax,)
        Array containing the sequence of approximations.
    """
    
    # Initialise the array with zeros
    p_array = np.zeros(Nmax)
    
    # Continue here:...
    if (f(a)*f(b))<0:
        for i in range(0,Nmax):
            p=(a+b)/2
            if f(a)*f(p)>0:
                a=p
                b=b
            elif f(b)*f(p)>0:
                b=p
                a=a
            p_array[i]=p
            i=i+1
    else:
        print("No roots found in the given interval")
    return p_array


#Question 3    
def fixedpoint_iteration(f,c,p0,Nmax):
    
    """
    Fixed Point Iteration Method: Returns a numpy array of the 
    sequence of approximations obtained by the fixed point iteration method.
    
    Parameters
    ----------
    f : function
        Input function for which the zero is to be found.
    c : real number
        Used in function g(x)
    p0 : real number
        Initial approximation.
    Nmax : integer
        Number of iterations to be performed.
        
    Returns
    -------
    p_array : numpy.ndarray, shape (Nmax,)
        Array containing the sequence of approximations.
    """
    
    # Initialise the array with zeros
    p_array = np.zeros(Nmax)
    
    # Continue here:...
    g=lambda x:x-(c*f(x))
    for i in range(0,Nmax):
        p=g(p0)
        print(p)
        p_array[i]=p
        p0=p
    return p_array


#Quesstion 4
def newton_method(f,dfdx,p0,Nmax):
    
        
    """
    Newton Method: Returns a numpy array of the 
    sequence of approximations obtained by the newton method.
    
    Parameters
    ----------
    f : function
        Input function for which the zero is to be found.
    dfdx : derivative of function
    p0 : real number
        Initial approximation.
    Nmax : integer
        Number of iterations to be performed.
        
    Returns
    -------
    p_array : numpy.ndarray, shape (Nmax,)
        Array containing the sequence of approximations.
    """
    
    # Initialise the array with zeros
    p_array = np.zeros(Nmax)
    
    # Continue here:...
    for i in range(0,Nmax):
        p=p0-(f(p0)/dfdx(p0))
        p_array[i]=p
        p0=p
    return p_array


#Question 5
def plot_convergence(p_exact,f,dfdx,c,p0,p1,Nmax,fig):
    x=np.arange(0,Nmax)
    y=bisection(f,p0,p1,Nmax)
    y1=fixedpoint_iteration(f,c,p0,Nmax)
    y2=newton_method(f,dfdx,p0,Nmax)
    y3=secant_method(f,p0,p1,Nmax)
    plt.xlabel("Number of Iterations")
    plt.ylabel("Error for various methods")
    plt.semilogy(x,np.abs(p_exact-y),'o',label="Bisection Method")
    plt.semilogy(x,np.abs(p_exact-y1),'o',label="Fixedpoint Iteration Method")
    plt.semilogy(x,np.abs(p_exact-y2),'o',label="Newton Method")
    plt.semilogy(x,np.abs(p_exact-y3),'o',label="Secant Method")
    plt.legend(loc="lower left",bbox_to_anchor=(0.4,0.2))
    plt.show() 
    
    
#Question 6
def secant_method(f,p0,p1,Nmax):
   
    
   """
   Secant Method: Returns a numpy array of the 
   sequence of approximations obtained by the secant method.
   
   Parameters
   ----------
   f : function
       Input function for which the zero is to be found.
   p0 : real number 
        Initial Approximation
   p1 : real number
        Second Approximation.
   Nmax : integer
       Number of iterations to be performed.
       
   Returns
   -------
   p_array : numpy.ndarray, shape (Nmax,)
       Array containing the sequence of approximations.
   """
   
   # Initialise the array with zeros
   p_array = np.zeros(Nmax)
   
   # Continue here:...
   for i in range(0,Nmax):
       if np.abs(f(p1)-f(p0))<(np.e)**-14:
           p_array[i]=p1 
       else:
           p=p1-(f(p1)*((p1-p0)/(f(p1)-f(p0))))
           p_array[i]=p
           p0=p1
           p1=p
   return p_array

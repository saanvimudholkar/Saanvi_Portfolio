"""
MATH2019 CW4 polynomial_interpolation module
"""

### Load other modules ###
import numpy as np
import matplotlib.pyplot as plt
### No further imports should be required ###

### Comment out incomplete functions before testing ###

#%%
def richardson(f,x0,h,k):
    """

    Parameters
    ----------
    f : Function of x
    x0 : Given Single Point
    h : Positive Integer
        Width of interval
    k : Positive Integer
        Level

    Returns
    -------
    deriv_approx : A float
                   Approximation of First derivative of f(x0)

    """
    array=np.zeros((k,k)) #Zero matrix of size(k,k)
    for i in range(0,k): #looping from 0 to k
        array[0,i]=(f(x0+(h/(2**i)))-f(x0-(h/(2**i))))/(2*(h/(2**i))) #Formula to find N1 given multiple of h
    for j in range(1,k): #looping from 1 to k
        for m in range(0,k-1): #looping from 0 to k-1
            ak=1/(1-(4**(j))) #Formula for ak
            bk=-((4**(j))/(1-(4**(j)))) #Formula for bk
            array[j,m]=ak*array[j-1,m]+bk*array[j-1,m+1] #Formula for Nk
    deriv_approx=array[k-1,0] #Approximate value of derivative
    return deriv_approx #return value to function

#%%
def richardson_errors(f,f_deriv,x0,n,h_vals,k_max):
    fig = plt.figure() #This line is required (once) before any other plotting commands
    error_matrix=np.zeros((k_max,n)) #Zero matrix of shape(k_max,n)
    for i in range(1,k_max+1): #Looping from 1 to k_max+1
         for j in range(0,n): #looping from 0 to n
            error_matrix[i-1,j]=np.abs(f_deriv(x0)-richardson(f,x0,h_vals[j],i)) #Formula for error calculation
    x=h_vals #Values of x
    for r in range(1,k_max+1): #looping from 1 to k_max+1
        y=error_matrix[r-1] #Values of y(error)
        plt.loglog(x,y,label="N(i)") #plotting log x against log y
    plt.xlabel("{h(i)} from i=0 to n-1") #labelling x axis
    plt.ylabel("Errors E(k-1,i) for i=0 to n-1") #labelling y axis
    plt.title("Error Terms at Different Levels") #labelling the graph
    plt.legend() #inserting legend
    plt.show() #producing graph
    return error_matrix, fig
    """

    Parameters
    ----------
    f : Function of x
    f_deriv : Derivative of function of x
    x0 : Initial condition
    n : Positive Integer
        Number of widths
    h_vals : a list of shape(n)
             Width 
    k_max : Positive Integer
            Maximum value of k

    Returns
    -------
    error_matrix : an np.ndarray of shape(k_max,n)
                   Error Terms
    fig : a matplotlib.figure Figure
    
   (a)    f(x)=sin(10x)
          ValueError is raised as during plotting the graph, x and y have 
          same first dimensions but different shapes i.e.(20,) and (4,)
          error_matrix is produced with the following values:
          [[1.66666663e-08 7.13555401e-08 3.05496782e-07 1.30793323e-06]
           [1.77635684e-15 0.00000000e+00 0.00000000e+00 1.42108547e-14]
           [0.00000000e+00 0.00000000e+00 1.77635684e-15 3.55271368e-15]
           [1.77635684e-15 0.00000000e+00 0.00000000e+00 0.00000000e+00]]
    (b)     f(x)=-x^2 for x<=0
          ValueError is raised as during plotting the graph, x and y have 
          same first dimensions but different shapes i.e.(20,) and (4,)
          error_matrix is produced with the following values:
          [[0. 0. 0. 0.]
           [0. 0. 0. 0.]
           [0. 0. 0. 0.]
           [0. 0. 0. 0.]]
          As x0 is zero and there is no constant term, we get a zero matrix.
    
    """
 


#### Your submission should have no code after this point ####




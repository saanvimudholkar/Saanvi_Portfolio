import numpy as np
import matplotlib.pyplot as plt
import newton #(run newton_solver(...) as newton.newton_solver(...))

def theta_ode_solver(a,b,f,df,N,y0,theta):
    """

    Parameters
    ----------
    a : Integer
        Starting Point of Interval
    b : Integer
        Endpoint of Interval
    f : Function of (t,y)
    df : Partial Derivative of f in terms of y
    N : Positive Integer
        Number of Interations
    y0 : Initial Condition
    theta : Positive Integer between 0 and 1

    Returns
    -------
    t : np.ndarray of shape(N+1,)
        Time points
    y : np.ndarray of shape(N+1,)
        Approximations at all time points

    """
    h=(b-a)/N #width of interval
    array=[y0] #initialising array with element y0
    p=[a] #initialising array with element a
    for i in range(0,N): #looping from 0 to N
        t=a #giving t value of a
        g=lambda x: x-y0-h*((1-theta)*f((t+h),x)+theta*f(t,y0)) #Defining function g(x)
        dg=lambda x: np.ones(np.shape(x)) #defining function g'(x)
        y=newton.newton_solver(g,dg,y0,100,1.0e-12) #Calling newton's solver
        y1=y0+h*(((1-theta)*f(t+h,y))+(theta*f(t,y0))) #Using theta schemes formula
        y0=y1 #switching values
        a=a+h #incrementing a by h
        p.append(a) #appending array p with new a
        array.append(y0) #appending array with y0 
    t=np.array(p) #making t an array of list p
    y=np.array(array) #making y and array of list array
    return t, y

def runge_kutta(a,b,f,N,y0,m,method):  
    """

    Parameters
    ----------
    a : Integer
        Starting point of interval
    b : Integer
        Endpoint of interval
    f : m-vector function
    N : Positive Integer
        Number of time steps
    y0 : Initial condition
    m : positive integer
    method : Runge Kutta method used
             (Method 1: Forward Euler Method
              Method 2: Midpoint Method
              Method 3: Heun's 3rd Order Method)

    Returns
    -------
    t : numpy.ndarray of shape(N+1,)
        Time points
    y : numpy.ndarray of shape(m,N+1)
        Approximations at all time points

    """
    t=np.zeros(N+1,) #Zero Matrix of shape (N+1)
    y=np.zeros((m,N+1)) #Zero matrix of shape(m,N+1)
    y[:,0]=y0 #1st column is values of y0
    h=(b-a)/N #width of interval
    for i in range(0,N+1): #looping from 0 to N+1
        t[i]=a+(i*h) #values of t
    if method==1: #Forward Euler Method
        r=a
        for i in range(1,N+1):
            y[:,i]=y0+h*f(r,y0)
            y0=y[:,i]
            r=a+i*h
    elif method==2: #Midpoint Method
        r=a
        for i in range(1,N+1):
            y[:,i]=y0+h*(f((r+(h/2)),(y0+((h/2)*f(r,y0)))))
            y0=y[:,i]
            r=a+i*h
    elif method==3: #Heun's 3rd Order Method
        r=a
        for i in range(1,N+1):
            k=f((r+(h/3)),(y0+((h/3)*f(r,y0))))
            g=f((r+((2*h)/3)),(y0+(((2*h)/3)*k)))
            y[:,i]=y0+(h/4)*(f(r,y0)+3*g)            
            y0=y[:,i]
            r=a+i*h
    return t, y




#MATH2019 CW3 polynomial_interpolation module


### Load other modules ###
import numpy as np
import matplotlib.pyplot as plt
### No further imports should be required ###

### Comment out incomplete functions before testing ###

#%%
def lagrange_poly(p,xhat,n,x,tol):
    """
    Parameters
    ----------
    p : Number of Nodal Points
    xhat : A numpy.ndarray consisting of p+1 Nodal Points
    n : Number of Evaluation Points
    x : A numpy.ndarray consisting of n Evalution Points 
    tol : Maximum Error Value allowed
    
    Returns
    -------
    lagrange_matrix : A numpy.ndarray of shape (n+1,p) where ijth entry is Li(xj)
    error_flag : Flag set to 0 if nodal points are distinct

    """
    lagrange_matrix=np.ones((p+1,n)) #Setting lagrange_matrix to be an array of size (p+1,n) with all entries equal to 1
    error_flag=0 #Setting error_flag to 0
    for i in range(0,p+1): #Looping i from 0 to p+1
        for j in range(0,p+1): #Looping j from 0 to p+1
            if i!=j: #Checking the condition whether i is not equal to j
                if np.abs(xhat[i]-xhat[j])<tol: #Checking whether the values of xhat are distinct
                    error_flag=1 #Setting error flag to one 1 if xhat values are not distinct
                lagrange_matrix[i,:]=lagrange_matrix[i,:]*((x-xhat[j])/(xhat[i]-xhat[j])) #Formula calculating lagrange_matrix
    return lagrange_matrix,error_flag #Returning lagrange_matrix and error flag
#%%
def uniform_poly_interpolation(a,b,p,n,x,f,produce_fig):
    """
    Parameters
    ----------
    a : Positive Integer
        Starting point of the interval
    b : Positive Integer
        End point of the interval
    p : Positive Integer
        Order of the polynomial
    n : Postive Integer
        Number of points in x
    x : np.ndarray consisting on n evaluating points
    f : Function f(x)
    produce_fig : Boolean expression

    Returns
    -------
    interpolant : np.ndarray of shape(n,)
                  Polynomial Interpolant values using uniformly distributed points
    fig : Graph of f(x) values and Polynomial Interpolant values at the same given x

    """
    fig = plt.figure() ##This line is required before any other plot commands
    interpolant=np.zeros((n,)) #A zero matrix of shape (n, )
    xhat=np.linspace(a,b,p+1) #Array of p+1 numbers in interval [a,b]
    matrix=lagrange_poly(p,xhat,n,x,1.0e-10)[0] #Calling function lagrange_poly and assigning matrix the value of lagrange_matrix from the called function 
    for i in range(0,p+1): #Looping from 0 to p+1
        interpolant=interpolant+(matrix[i,:]*f(xhat[i])) #Calculating the polynomial interpolant
    if produce_fig==True: #Checking condition if the Boolean Expression produce_fig is True 
        y1=f(x) #Assigning y1 with all values of f(x) for all given x
        y2=interpolant #Assigning y2 with the values of polynomial interpolant calculated
        plt.plot(x,y1, label="f(x)") #Plotting x VS f(x)
        plt.plot(x,y2, label="Interpolant Values") #Plotting x VS Interpolant values 
        plt.xlabel("Values of x") #Labelling the x axis
        plt.legend() #Declaring a legend
        plt.title("Plot showing the actual values of f(x) and interpolant values") #Giving a title to the graph
    else:
        return interpolant, None
    return interpolant, fig #Returning values of interpolant and plotting the figure

#%%
def nonuniform_poly_interpolation(a,b,p,n,x,f,produce_fig):
    """
    Parameters
    ----------
    a : Positive Integer
        Starting point of interval        
    b : Posiitve Integer
        End point of interval
    p : Positive Integer
        Order of polynomial
    n : Positive Integer
        Number of points in x
    x : np.ndarray consisting of n evaluating points
    f : Function f(x)
    produce_fig : Boolean expression for producing graph

    Returns
    -------
    nu_interpolant : np.ndarray of shape(n,)
                     Polynomial Interpolant values calculated using non-uniform points
    fig : Graph of f(x) values and Polynomial Interpolant values at the same given x

    """
    fig=plt.figure() #Declaring a figure
    xhat=np.zeros(p+1,) #Zero array of shape (p+1)
    y=np.zeros(p+1,) #Zero array of shape(p+1)
    nu_interpolant=np.zeros((n,)) #Zero array of shape(n)
    for i in range(0,p+1): #Looping from 0 to p+1
        xhat[i]=np.cos((((2*i)+1)/(2*(p+1)))*np.pi) #Calculating xhat values at given i values
    m=((b-a)/2) #Width of each interval
    for k in range(0,p+1): #Looping from 0 to p+1
        y[k]=m*(xhat[k]+1)+a #Linearly transating obtained values of xhat
    matrix=lagrange_poly(p,y,n,x,1.0e-10)[0] #Calling the function lagrange_poly
    for j in range(0,p+1): #Looping from 0 to p+1
        nu_interpolant=nu_interpolant+(matrix[j,:]*f(y[j])) #Calculating polynomial interpolating values for xhat values calculated
    if produce_fig==True: #Checking condition of produce_fig
        y1=f(x) #assigning f(x) to y1
        y2=nu_interpolant #assigning nu_interpolant to y2
        plt.plot(x,y1,label="f(x)") #Plot x vs y1
        plt.plot(x,y2,label="nu-interpolant values") #Plot x vs y2
        plt.xlabel("Values of x") #Labelling x axis
        plt.title("Graph of f(x) and Interpolant values for the same given x") #Labelling the graph
        plt.legend() #Inserting legend
    return nu_interpolant, fig #returning values and figure

#%%
def compute_errors(a,b,n,P,f):
    """
    Parameters
    ----------
    a : Positive Integer
        Starting point of the interval
    b : Positive Integer
        End point of the interval
    n : Positive integer
        Number of evaluating points
    P : np.ndarray of shape n
        Array of degrees of polynomial
    f : Function f(x)


    Returns
    -------
    error_matrix : np.ndarray of shape(2,n)
                   Array of error terms for uniform polynomial interpolation and non-uniform polynomial interpolation
    fig : Graph of x vs Error terms of uniform points and x vs error terms of non-uniform points

    """
    fig=plt.figure() #Declaring a figure
    x=np.linspace(a,b,2000) #Array of 2000 uniformly spaced points in the given interval [a,b]
    array1=[] #Declaring array1
    array2=[] #Declaring array2
    for i in range(0,len(P)): #Looping from 0 to the lenght of P
        matrix1=uniform_poly_interpolation(a,b,P[i],2000,x,f,False)[0] #Calling function uniform_poly_interpolation
        matrix2=nonuniform_poly_interpolation(a,b,P[i],2000,x,f,False)[0] #Calling function nonuniform_poly_interpolation
        max1=max(np.abs(matrix1-f(x))) #Finding maximum value of error term in uniform sequence
        max2=max(np.abs(matrix2-f(x))) #Finding maximum value of error term in non uniform sequence
        array1.append(max1) #Appending array1 with value of max1
        array2.append(max2) #Appending array2 with value of max2
    b=[array1,array2] #Concatening array1 and array2
    error_matrix=np.asarray(b) #Converting b into an array
    plt.semilogy(P,error_matrix[0],label="Error terms using uniform polynomial interpolation") #Plotting Degrees of polynomials vs errors in uniform interpolation 
    plt.semilogy(P,error_matrix[1],label="Error terms using non-uniform polynomial interpolation") #Plotting Degrees of polynomials vs errors in nonuniform interpolation
    plt.title("Error Terms") #Giving a title to the graph
    plt.legend() #Showing legend
    return error_matrix,fig #returning error_matrix and graph
#%%
def piecewise_interpolation(a,b,p,m,n,x,f,produce_fig):
    """
    Parameters
    ----------
    a : Postive Integer
        Starting of the interval
    b : Positive Integer
        End of the interval
    p : Positive Integer
        Order of Polynomial
    m : Positive integer
        Number of subintervals
    n : Positive integer
        Number of points
    x : np.ndarray of shape n
        Evaluating points
    f : Function f(x)
    produce_fig : Boolean expression for graph

    Returns
    -------
    p_interpolant : np.ndarray
                    piecewise interpolant values
    fig : Graph of x vs f(x) and x vs Piecewise Interpolant values

    """
    p_interpolant=[]
    fig=plt.figure()
    div=np.linspace(a,b,m+1)
    a=[]
    b=[]
    for i in range(0,m):
        arr=nonuniform_poly_interpolation(div[i],div[i+1],p,n,x,f,False)[0]
        a.append(arr)  
    for j in range(0,m):
        for k in range(0,n):
            if j<k<=(j+1):
                b.append(a[j][k])
    b.insert(0,a[0][0])
    b.insert(n-1,a[m-1][n-1])
    p_interpolant=b
    if produce_fig==True:
        plt.plot(x,p_interpolant)
        plt.plot(x,f(x))
    else:
        return p_interpolant,None
    return p_interpolant,fig

#### Your submission should have no code after this point ####


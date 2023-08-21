"""
MATH2050 CW2 systemsolvers module

@author: Saanvi Mudholkar
"""

import numpy as np
import matplotlib.pyplot as plt
import backward as bw

def no_pivoting(A,b,n,c):
    
    """
    Returns the augmented matrix M arrived at by starting from the augmented
    matrix [A b] and performing forward elimination without row interchanges
    until all of the entries below the main diagonal in the first c columns
    are 0.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the matrix A in the linear system Ax=b.
    b : numpy.ndarray of shape (n,1)
        array representing the vector b in the linear system Ax=b.
    n : integer
        positive integer.
    c : integer
        positive integer that is at most n-1.
    
    Returns
    -------
    M : numpy.ndarray of shape (n,n+1)
        2-D array representing the matrix M.
    """
    
    # Create the initial augmented matrix
    M=np.hstack((A,b))
    
    # Continue here:...
    for i in range(c):
        if M[i,i] == 0.0:
            print('Division by zero')
        for j in range(i+1, n):
            ratio = M[j,i]/M[i,i]
            for k in range(n+1):
                M[j,k]= M[j,k] - ratio * M[i,k] 
    return M

def no_pivoting_solve(A,b,n):
    
    """
    Returns the solution x to Ax = b computed using Gaussian elimination without 
    pivoting.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the matrix A in the linear system Ax=b.
    b : numpy.ndarray of shape (n,1)
        array representing the vector b in the linear system Ax=b.
    n : integer
        positive integer.
    
    Returns
    -------
    X : numpy.ndarray of shape (n,1)
        2-D array representing the solution to matrix M.
    """
    M=no_pivoting(A,b,n,n-1)
    X=bw.backward_substitution(M, n)
    return X

#Question 2
                      
def find_max (M,n,i ):
    """
    Returns the row index for maximum value in the given row.
    
    Parameters
    ----------
    M : numpy.ndarray of shape (n,n+1)
        array representing the augumented matrix M.
    n : integer
        positive integer.
    i: given row index
    Returns
    -------
    r : row index for maximum value in given row.
    """
    max=0
    r=0
    for j in range(0,n):
        if M[i,j]>max:
           max=M[i,j]
           r=j
    return r

def partial_pivoting(A,b,n,c):
    """
    Returns the solution after performing partial pivoting by row swapping.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the matrix A in the linear system Ax=b.
    b : numpy.ndarray of shape (n,1)
        array representing the vector b in the linear system Ax=b.
    n : integer
        positive integer.
    c : integer
        positive integer that is at most n-1.
    
    Returns
    -------
    M : numpy.ndarray of shape (n,n+1)
        2-D array representing the matrix M after partial pivoting.
    """
    M=np.hstack((A,b))
    r=find_max(M.T,n,0)
    t=[]
    for p in range(0,n+1):
        t=M[0,p]
        M[0,p]=M[r,p]
        M[r,p]=t
    for s in range(c):
        if M[s,s] == 0.0:
            print('Division by zero')
        for j in range(s+1, n):
            ratio = M[j,s]/M[s,s]
            for k in range(n+1):
                M[j,k]= M[j,k] - ratio * M[s,k]
    return M

def partial_pivoting_solve(A,b,n):
    M=np.hstack((A,b))
    M=partial_pivoting(A,b,n,n-1)
    X=bw.backward_substitution(M,n)
    return X

#Question 3

def Doolittle(A,n):
    """
    Returns the Lower Triangular and Upper Triangular matrix in the 
    Dolittle factorisation.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the matrix A in the linear system Ax=b.
    n : integer
        positive integer.
        
    Returns
    -------
    L : numpy.ndarray of shape (n,n)
        2-D array representing the Lower Triangular Matrix for 
        Dolittle Factorisation.
    U : numpy.ndarray of shape (n,n)
        2-D array representing the Upper Triangular Matrix for 
        Dolittle Factorisation.
    """
    L=np.identity(3)
    U=A
    for j in range(n-1):
        for i in range(j+1,n):
            L[i,j]=(A[i,j]-sum(U[k,j]*L[i,k] for k in range(1,j-1)))/U[j,j]
    for i in range(n-1):
        if U[i,i] == 0.0:
            print('Division by zero')
        for j in range(i+1, n):
            ratio = U[j,i]/U[i,i]
            for k in range(n):
                U[j,k]= U[j,k] - ratio * U[i,k] 
    return L,U   

#Question 4

def Gauss_Seidel(A,b,n,x0,tol,maxits):
    """
    Returns the approximations of solved linear system using Gauss-Seidel 
    method.
    
    Parameters
    ----------
    A : numpy.ndarray of shape (n,n)
        array representing the matrix A in the linear system Ax=b.
    b : numpy.ndarray of shape (n,1)
        array representing the matrix b in the linear system Ax=b.
    n : integer
        positive integer.
   x0 : numpy.ndarray of shape (n,1)
        array representing initial approximation.
  tol : positive real number.
maxits: positive integer
                
    Returns
    -------
    p : numpy.ndarray of shape (n,1)
         array represnting approximation of system after maxits iterations.
    """
    D=np.zeros(shape=(3,3))
    L=np.zeros(shape=(3,3))
    U=np.zeros(shape=(3,3))
    for i in range(0,n):
        D[i,i]=A[i,i]
    for i in range(0,n-1):
        for j in range(i+1,n):
            L[j,i]=np.abs(A[j,i])
            U[i,j]=np.abs(A[i,j])
    R=D-L
    I=np.linalg.inv(R)
    T=np.dot(I,U)
    c=np.dot(I,b)
    for i in range(1,maxits+1):
        xi=(np.dot(T,x0))+c
        x0=xi
        p=xi
    if i<=maxits and np.max(np.abs(b-np.matmul(A,x0)))<tol :
        return p
    else:
        p="Desired tolerance not reached after maxits iterations have been performed"
        return p
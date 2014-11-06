import numpy as np

"""
A2-Part-3: Implement the discrete Fourier transform (DFT)

Write a function that implements the discrete Fourier transform (DFT). Given a sequence x of length
N, the function should return its DFT, its spectrum of length N with the frequency indexes ranging from 0 
to N-1.

The input argument to the function is a numpy array x and the function should return a numpy array X which 
is of the DFT of x.

EXAMPLE: If you run your function using x = np.array([1, 2, 3, 4]), the function shoulds return the following numpy array:
array([10.0 + 0.0j,  -2. +2.0j,  -2.0 - 9.79717439e-16j, -2.0 - 2.0j])

Note that you might not get an exact 0 in the output because of the small numerical errors due to the
limited precision of the data in your computer. Usually these errors are of the order 1e-15 depending
on your machine.
"""
def DFT(x):
    """
    Input:
        x (numpy array) = input sequence of length N
    Output:
        The function should return a numpy array of length N
        X (numpy array) = The N point DFT of the input sequence x
    """
    ## Your code here

    N = x.size
    points = np.empty((N, N), dtype=np.complex128)
    for n in xrange(N):
        points[n] = n * np.arange(N)
    arg = 2 * np.pi * points / N
    bases = np.cos(arg) - 1j * np.sin(arg)
    return np.dot(bases, x)

def DFT2(x):
    """
    Input:
        x (numpy array) = input sequence of length N
    Output:
        The function should return a numpy array of length N
        X (numpy array) = The N point DFT of the input sequence x
    """
    ## Your code here

    # This version is slower, it doesn't take full advantage of
    # numpy's vectorized operations

    N = x.size
    arr = np.empty(N, dtype=np.complex128)
    print x
    for k in xrange(N):
        basis = genComplexSine(k, N)
        arr[k] = np.dot(x, basis)
    return arr

def get_points(N):
    points = np.empty((N, N))
    for n in xrange(N):
        points[n] = (n+1) * np.arange(N)
    return points


def genComplexSine(k, N):
    """
    Inputs:
        k (integer) = frequency index of the complex sinusoid of the DFT
        N (integer) = length of complex sinusoid in samples
    Output:
        The function should return a numpy array
        cSine (numpy array) = The generated complex sinusoid (length N)
    """
    ## Your code here
    samples = np.arange(N)
    arg = 2 * np.pi * k * samples / N
    print arg
    return np.cos(arg) - 1j * np.sin(arg)

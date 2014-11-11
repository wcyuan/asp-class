from __future__ import division

import numpy as np
from scipy.signal import get_window
import math
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import dftModel as DFT
import utilFunctions as UF
import cPickle

""" 
A5-Part-1: Minimizing the frequency estimation error of a sinusoid

Write a function that estimates the frequency of a sinusoidal signal at a given time instant. The 
function should return the estimated frequency in Hz, together with the window size and the FFT 
size used in the analysis.  

The input arguments to the function are the wav file name including the path (inputFile) containing 
the sinusoidal signal, and the frequency of the sinusoid in Hz (f). The frequency of the input sinusoid  
can range between 100Hz and 2000Hz. The function should return a three element tuple of the estimated 
frequency of the sinusoid (fEst), the window size (M) and the FFT size (N) used.

The input wav file is a stationary audio signal consisting of a single sinusoid of length 1 second. 
Since the signal is stationary you can just perform the analysis in a single frame, for example in 
the middle of the sound file (time equal to .5 seconds). The analysis process would be to first select 
a fragment of the signal equal to the window size, M, centered at .5 seconds, then compute the DFT 
using the dftAnal function, and finally use the peakDetection and peakInterp functions to obtain the 
frequency value of the sinusoid.

Use a Blackman window for analysis and a magnitude threshold t = -40 dB for peak picking. The window
size and FFT size should be chosen such that the difference between the true frequency (f) and the 
estimated frequency (fEst) is less than 0.05 Hz for the entire allowed frequency range of the input 
sinusoid. The window size should be the minimum positive integer of the form 100*k + 1 (where k is a 
positive integer) for which the frequency estimation error is < 0.05 Hz. For a window size M, take the
FFT size (N) to be the smallest power of 2 larger than M. 

HINT: If the specified frequency range would have been 440-8000 Hz, the parameter values that satisfy 
the required conditions would be M = 1101, N = 2048. Note that for a different frequency range, like 
the one specified in the question, this value of M and N might not work. 

"""

DEFAULT_WINDOW = 'blackman'
DEFAULT_FS = 44100
DEFAULT_THRESHOLD = -40 # decibel
DEFAULT_TIME = 1 # second
DEFAULT_LOW_FREQ = 100
DEFAULT_HIGH_FREQ = 2000
DEFAULT_FREQ_ERROR = 0.05

def minFreqEstErr(inputFile, f):
    """
    Inputs:
            inputFile (string) = wav file including the path
            f (float) = frequency of the sinusoid present in the input audio signal (Hz)
    Output:
            fEst (float) = Estimated frequency of the sinusoid (Hz)
            M (int) = Window size
            N (int) = FFT size
    """
    print "freq: {0}, inputFile {1}".format(f, inputFile)

    # analysis parameters:
    window = 'blackman'
    t = -40
    
    ### Your code here
    (fs, x) = UF.wavread(inputFile)

    M = 2101
    (mX, pX, ploc, iploc, ipmag, ipphase, fEst, N) = run_one_estimate(
        x, fs, M, window=window, t=t)

    print fEst
    return (fEst[0], M, N)

def min_power_2(M):
    """
    Return the first power of 2 higher than M
    """
    N = 1
    while N < M:
        N *= 2
    return N

def run_one_estimate(x, fs, M, window=DEFAULT_WINDOW, t=DEFAULT_THRESHOLD):
    center_sample = int(len(x) / 2)
    start_sample = center_sample - int(M/2)
    end_sample = start_sample + M
    N = min_power_2(M)
    x1 = x[start_sample:end_sample]
    w = get_window(window, M)
    mX, pX = DFT.dftAnal(x1, w, N)
    ploc = UF.peakDetection(mX, t)
    iploc, ipmag, ipphase = UF.peakInterp(mX, pX, ploc)
    fEst = iploc * fs / N
    return (mX, pX, ploc, iploc, ipmag, ipphase, fEst, N)

def find_best_window_size(x, f, fs, window=DEFAULT_WINDOW, t=DEFAULT_THRESHOLD, err=DEFAULT_FREQ_ERROR, verbose=False):
    """
    Find the first window size of the form (100K+1) that is within the
    error for a particular frequency.

    Unused.  Thought this was what we were supposed to do, but it isn't.
    """
    for k in xrange(1, 1000):
        M = 100 * k + 1
        (mX, pX, ploc, iploc, ipmag, ipphase, fEst, N) = run_one_estimate(
            x, fs, M, window=window, t=t)
        if verbose:
            print "f={0} M={1} N={2} fEst={3} error={4}".format(f, M, N, fEst, f-fEst)
        if iploc.size > 0 and all(abs(fEst - f) < err):
            return (fEst, M, N)

def genSine(A, f, phi, fs, t):
    """
    Inputs:
        A (float) =  amplitude of the sinusoid
        f (float) = frequency of the sinusoid in Hz
        phi (float) = initial phase of the sinusoid in radians
        fs (float) = sampling frequency of the sinusoid in Hz
        t (float) =  duration of the sinusoid (is second)
    Output:
        The function should return a numpy array
        x (numpy array) = The generated sinusoid
    """
    ## Your code here
    
    return A * np.cos(2 * np.pi * f * np.arange(fs * t) / fs + phi)

def find_best_window_over_range(low_f=DEFAULT_LOW_FREQ, high_f=DEFAULT_HIGH_FREQ, fs=DEFAULT_FS,
                                err=DEFAULT_FREQ_ERROR, t=DEFAULT_THRESHOLD, window=DEFAULT_WINDOW,
                                verbose=False):
    """
    For each frequency, find the first window_size of the form 100*K+1
    that is within the error.  Then find the highest of all those window
    sizes.

    Unused.  Thought this was what we were supposed to do, but it isn't.
    """
    best_M = 0
    for f in xrange(low_f, high_f+1):
        x = genSine(A=1, f=f, phi=0, fs=fs, t=1)
        (fEst, M, N) = find_best_window_size(x, f, fs, window=window, t=t, err=err, verbose=verbose)
        print "Freq {0} requires M={1} N={2} fEst={3}".format(f, M, N, fEst)
        if M > best_M:
            best_M = M
            print "Best window size is now {0}".format(M)
    return best_M


def find_best_window_for_all_freq(low_f=DEFAULT_LOW_FREQ, high_f=DEFAULT_HIGH_FREQ, fs=DEFAULT_FS,
                                  err=DEFAULT_FREQ_ERROR, t=DEFAULT_THRESHOLD, window=DEFAULT_WINDOW,
                                  verbose=True):
    """Find the lowest window size of the form 100K+1 that allows all
    frequencies in the range to have error less than the given error

    It's not very clear from the instructions, but this is actually
    what you're supposed to run.  Run this once to figure out the
    window size to use, then hard code that window size in the main
    function.
    """
    best_M = 0
    for k in xrange(1, 1000):
        M = 100 * k + 1
        for f in xrange(low_f, high_f+1):
            x = genSine(A=1, f=f, phi=0, fs=fs, t=1)
            (mX, pX, ploc, iploc, ipmag, ipphase, fEst, N) = run_one_estimate(
                x, fs, M, window=window, t=t)
            if verbose:
                print "f={0} M={1} N={2} fEst={3} error={4}".format(f, M, N, fEst, f-fEst)
            if iploc.size <= 0 or any(abs(fEst - f) > err):
                break
        else:
            return (M, N)



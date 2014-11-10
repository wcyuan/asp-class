from __future__ import division

import os
import sys
import numpy as np
from scipy.signal import get_window
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import stft
import utilFunctions as UF


eps = np.finfo(float).eps

"""
A4-Part-3: Computing band-wise energy envelopes of a signal

Write a function that computes band-wise energy envelopes of a given audio signal by using the STFT.
Consider two frequency bands for this question, low and high. The low frequency band is the set of all the 
frequencies from 0 - 3000 Hz and the high frequency band is the set of all the frequencies from 
3000 - 10000 Hz. At a given frame, the value of the energy envelope of a band can be computed as the 
sum of squared values of all the frequency coefficients in that band. Compute the energy envelopes in 
decibels. 

Refer to "A4-STFT.pdf" document for further details on computing bandwise energy.

The input arguments to the function are the wav file name including the path (inputFile), window 
type (window), window length (M), FFT size (N) and hop size (H). The function should return a numpy 
array with two columns, where the first column is the energy envelope of the low frequency band and 
the second column is that of the high frequency band.

Use stft.stftAnal() to obtain the STFT magnitude spectrum for all the audio frames. Then compute two 
energy values for each frequency band specified. While calculating frequency bins for each frequency 
band, consider only the bins that are within the specified frequency range. For example, for the low 
frequency band consider only the bins with frequency > 0 Hz and < 3000 Hz. This way we also remove the 
DC offset in the signal in energy envelope computation.

To get a better understanding of the energy envelope and its characteristics you can plot the envelopes 
together with the spectrogram of the signal. You can use matplotlib plotting library for this purpose. 
To visualize the spectrogram of a signal, a good option is to use colormesh. You can reuse the code in
sms-tools/lectures/4-STFT/plots-code/spectrogram.py. Either overlay the envelopes on the spectrogram 
or plot them in a different subplot. Make sure you use the same range of the x-axis for both the 
spectrogram and the energy envelopes.

EXAMPLE: Running your code on piano.wav file with window = 'blackman', M = 513, N = 1024, H = 128, in 
the plots you can clearly notice the sharp attacks and decay of the piano notes (See figure in the 
accompanying pdf). In addition, you can also visually analyse which of the two energy envelopes is better 
for detecting onsets of the piano notes.
"""
def computeEngEnv(inputFile, window, M, N, H):
    """
    Inputs:
            inputFile (string): input sound file (monophonic with sampling rate of 44100)
            window (string): analysis window type (choice of rectangular, triangular, hanning, hamming, 
                blackman, blackmanharris)
            M (integer): analysis window size (odd positive integer)
            N (integer): FFT size (power of 2, such that N > M)
            H (integer): hop size for the stft computation
    Output:
            The function should return a numpy array engEnv with shape Kx2, K = Number of frames
            containing energy envelop of the signal in decibles (dB) scale
            engEnv[:,0]: Energy envelope in band 0 < f < 3000 Hz (in dB)
            engEnv[:,1]: Energy envelope in band 3000 < f < 10000 Hz (in dB)
    """
    
    ### your code here

    fs, x = UF.wavread(inputFile)
    w = get_window(window, M)
    (mX, pX) = stft.stftAnal(x, fs, w, N, H)

    numFrames = int(mX[:,0].size)
    frmTime = H*np.arange(numFrames)/float(fs)
    binFreq = np.arange(N/2+1)*float(fs)/N

    cutoff1 = 3000
    cutoff2 = 10000

    cutoff_bucket1 = np.ceil(float(cutoff1) * N / fs)
    cutoff_bucket2 = np.ceil(float(cutoff2) * N / fs)

    low_band = mX[:,1:cutoff_bucket1]
    high_band = mX[:,cutoff_bucket1:cutoff_bucket2]

    E = np.zeros((numFrames, 2))
    E[:,0] = by_frame_energy(low_band)
    E[:,1] = by_frame_energy(high_band)


    #plot_energies(mX, fs, inputFile, M, N, H, E)

    return E

def by_frame_energy(band):
    """
    @param band: a K x F numpy array which is the output of an STFT where K
    is the number of frames and F is the number of frequency buckets
    """
    # convert from decibels to linear
    linear = np.power(10, band/20.0)
    # compute energy (sum of squares along the frequency buckets, i.e. along the second axis)
    squares = linear * linear
    energy = np.sum(squares, axis=1)
    # return to decibels
    return 10*np.log10(energy)

def plot_energies(mX, fs, inputFile, M, N, H, E):
    plt.figure(1, figsize=(9.5, 6))
    
    plt.subplot(211)
    numFrames = int(mX[:,0].size)
    frmTime = H*np.arange(numFrames)/float(fs)
    binFreq = np.arange(N/2+1)*float(fs)/N
    plt.pcolormesh(frmTime, binFreq, np.transpose(mX))
    plt.title('mX ({3}), M={0}, N={1}, H={2}'.format(M, N, H, inputFile))
    plt.autoscale(tight=True)
    
    plt.subplot(212)
    plt.plot(frmTime, E[:,0])
    plt.plot(frmTime, E[:,1])
    plt.autoscale(tight=True)
    
    plt.tight_layout()
    plt.legend(("low", "high"))
    #plt.savefig('spectrogram.png')
    #plt.show()



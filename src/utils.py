"""
Useful functions for the device operation, this file is based
in the core code developed for the acoustic monitoring station.
"""

import os
import numpy as np
from scipy.signal import (
    lfilter, sosfilt, butter, bilinear
    )
from scipy.io import wavfile

# Current file directory
DIR = os.path.dirname(os.path.realpath(__file__))

## Constant definitions:
# - - - - - - - - - - - - - - - 

# sampling frequency [Hz]
FS = 48000
# reference frequency [Hz]
FR = 1000

## A and C weighting:
## zero-pole-gain system from IEC 61672

# C weighting
Z_C = np.array([0, 0])
P_C = np.array([-2*np.pi*20.598997057568145, 
                -2*np.pi*20.598997057568145,
                -2*np.pi*12194.21714799801,
                -2*np.pi*12194.21714799801])
K_C = (10**(0.062/20))*P_C[3]**2

# A weighting
Z_A = np.append(Z_C, [0, 0])
P_A = np.insert(P_C, 
                [2, 2], 
                [-2*np.pi*107.65264864304628, 
                -2*np.pi*737.8622307362899])
K_A = (10**(2/20))*P_A[4]**2

## Loading coefficients:
# - - - - - - - - - - - - - - - 

H_inv = np.load(DIR + '\H_inv.npy')
sos_A = np.load(DIR + '\sos_A.npy')
sos_C = np.load(DIR + '\sos_C.npy')

## Functions:
# - - - - - - - - - - - - - - -

def wavread(x):
    """ Reads and normalizes the input WAV file 'x' based on its location. """
    norm_fact = {'int16': (2**15)-1, 'int32': (2**31)-1, 
                 'int64': (2**63)-1, 'float32': 1.0, 
                 'float64': 1.0}
    fs, y = wavfile.read(x)
    y = np.float32(y)/norm_fact[y.dtype.name]
    return fs, y

def spec(x, nfft=2**12, fs=FS):
    """ Compute the FFT spectrum of the input audio signal 'x'. 
    The number of points for which to compute the FFT can be 
    specified by 'nfft'. """
    freq = fs/2*np.linspace(0, 1, int(nfft/2)+1)  # FFT frequency vector
    sp = np.fft.rfft(x, nfft)
    spp = np.abs(sp)/nfft
    return spp, freq

def ref_f(X, f):
    """ References the input signal 'X' at 1 kHz according
    to the frequency array 'f'. """
    idx = (np.abs(f-1000)).argmin()
    return X-X[idx]

def fft_average(x, nfft=2**16):
    """ Receives a tuple or list of audio arrays 'x', calculates and 
    normalizes the FFT spectrum of each of them and returns their average. 
    With 'nfft' the number of FFT points can be modified. """
    X = np.empty([len(x), int(nfft/2+1)])
    xn = np.empty([len(x), int(nfft/2+1)])
    y = np.empty(int(nfft/2))
    for i in range(len(x)):
        X[i,:], freq = spec(x[i], nfft=nfft)
        xn[i,:] = 1 + ref_f(X[i,:], freq)
    y = np.apply_along_axis(np.mean, 0, xn)
    return y, freq

def smooth_spectrum(X, freq, Noct=3):
    """ Apply a 1/Noct smoothing to the input frequency spectrum 'X'
    at the 'freq' frequencies.
    Calculates the i-th smoothed spectral coefficient 'X_oct(i)' as 
    the sum of the windowed spectrum. The window is a Gaussian whose
    center frequency is 'freq(i)', and whose sigma is proportional
    to 'freq(i)/Noct'. """
    # initial spectrum
    X_oct = X.copy()
    if Noct > 0:
        for i in np.arange(np.where(freq>0)[0][0], len(freq)):
            g = gauss_f(freq=freq, f0=freq[i], Noct=Noct)
            X_oct[i] = sum(g*X)
    # avoid undershoot
    if min(X) >= 0:
        X_oct[X_oct<0] = 0
    return X_oct

def gauss_f(freq, f0, Noct):
    """ Calculate frequency-domain Gaussian with unity gain. 
    Inputs:
    freq: array of frequencies in which to do the computation.
    f_0: center frequency.
    Noct: determines the bandwith as f_0/Noct. """
    # standard deviation
    sigma = (f0/Noct)/np.pi    
    # Gaussian
    g = np.exp(-(((freq-f0)**2)/(2*(sigma**2))))
    # magnitude normalizaition
    return g/sum(g)

def inverse_filter(x):
    """ Appyling the system's inverse filter, according to the
    coefficients from the designed included on Notebook 1. """
    return lfilter(H_inv, 1, x)

def bilinear_zpk(z, p, k, fs=FS):
    """ Returns the zero-pole-gain system in the z-domain
    from the same system in the s-domain through the bilinear
    transform. """
    deg = len(p) - len(z)
    fs2 = 2.0*fs
    # Bilinear transform
    z_b = (fs2 + z)/(fs2 - z)
    p_b = (fs2 + p)/(fs2 - p)
    z_b = np.append(z_b, -np.ones(deg))
    k_b = k*np.real(np.prod(fs2 - z)/np.prod(fs2 - p))
    return z_b, p_b, k_b

def filt_A(x):
    """ Applies the A weighting filter to the input signal 'x' """
    return sosfilt(sos_A, x)

def filt_C(x):
    """ Applies the C weighting filter to the input signal 'x' """
    return sosfilt(sos_C, x)

def upp_low_freqs(f_mid, Noct=1.0):
    """ Returns the upper and lower Noct upper and lower band-edge
    frequencies according to the given 'f_mid' mid-band frequencies."""
    return np.around(f_mid*Noct**(-1/(2*Noct)), 5), np.around(f_mid*Noct**(1/(2*Noct)), 5)

def sos_bp(low, upp, fs=FS, order=4):
    """ Returns the coefficients to a second-order-section 
    Butterworth band-pass filter. The filter is built with 
    'order' orders and according to the upper and lower 
    edge-band frequencies ('low' and 'up'). """
    nyq = 0.5*fs
    f1 = low/nyq
    f2 = upp/nyq
    sos = butter(order, [f1, f2], btype='band', output='sos')
    return sos

def bp_filt(x, low, upp, fs=FS, order=4):
    """ Applies a second-order-section band-pass filter 
    to the input signal 'x'.
    'low' and 'upp', the upper and lower edge-band frequencies
    must be specified. """
    sos = sos_bp(low, upp, fs=fs, order=order)
    return sosfilt(sos, x)

def generate_weighting_filter(t=1.0, fs=FS):
    """ Generates the coefficients for a second-order-section
    'Slow' or 'Fast' time-weighting filter, according to the
    time constant 't' and the sampling frequency 'fs'. """
    return butter(1, 1/(np.pi*t*fs), output='sos')

def time_weighting(x, which='Slow', fs=FS):
    """ Applies time-weigthing to the input signal 'x'
    according to the sampling frequecy 'fs'. 
    The argument 'which' detemines if 'Slow' (default)
    or 'Fast' time constants will be applied. """
    # defining time constant
    t = 1.0
    if which == 'Fast':
        t = 0.125
    # computing the coefficients
    sos = generate_weighting_filter(t=t, fs=fs)
    # appliying the filter
    return np.sqrt(sosfilt(sos, x**2))

def rms(x):
    """ RMS level calculatiaon for the input signal 'x'. """
    return np.sqrt(np.mean(x)**2)

def rms_t(x, t=1.0, fs=FS):
    """ RMS level calculation for the input signal 'x',
    in fragments of length 't'.
    It's expected to recieve a 1-dimension 'x' signal
    or a 2-dimension 'x' singal with x.shape[0] beign the 
    amount of frequency bands which the signal is comprised of. 
    The output will have the same number of dimensions as the input. """
    x = np.array(x, dtype='float64')
    N = int(np.floor(t*fs))
    if x.ndim == 1 :
        start = np.arange(0, len(x), N) 
        end = np.append(start[1:], x.shape[0])
        y= np.empty(len(start))
        for i in range(len(start)):
            y[i] = rms(x[start[i]:end[i]])
    else:
        start = np.arange(0, x.shape[0], N) 
        end = np.append(start[1:], x.shape[0])
        y= np.empty([x.shape[0], len(start)])
        for i in range(x.shape[0]):
            for j in range(len(start)):
                y[i, j] = rms(x[i, start[j]:end[j]])
    return y

def mean_db(x):
    """ Get the energetic average of the input array of sound
    levels in dB'X'. """
    return 10*np.log10(np.mean(10**(x/10)))
"""
Useful functions for the device operation, this file is based
in the core code developed for the acoustic monitoring station
and will be explained and used throughout the notebooks.

In this file we will include:

    - 
"""


def spec(x, nfft=2**16, window='rect'):
    """ Compute the FFT spectrum of the input audio signal 'x'. 
    The number of points for which to compute the FFT can be 
    specified by 'nfft'. Its possible to select a Hanning or Flattop
    window, in which case corrections are applied to obtain accurate 
    measurements based on the spectrum. """

    freq = fs/2*np.linspace(0, 1, int(nfft/2)+1)  # FFT frequency vector
   
    if window == 'hann':
        window = np.hanning(len(x))  # Hann window
        CPG = 6.02 # Coherent Power Gain [dB]
        data = x*window

    elif window == 'flat':
        window = flattop(len(x))  # flattop window
        CPG = 13.3 # Coherent Power Gain [dB]
        data = x*window

    else:
        data = x

    sp = np.fft.rfft(data, nfft)
    spp = np.abs(sp)/nfft

    return spp, freq


def ref_f(x, f):
    """ References the input signal 'x' at 1 kHz according
    to the frequency array 'f'. """

    idx = (np.abs(f-1000)).argmin()

    return x-x[idx]


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
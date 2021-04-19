# My engineering thesis

In this repo you'll find some Jupyter notebook examples of different funcionality that was used for my **engineering thesis**: **_"Development of a low-cost urban acoustic monitoring station"_**.

The core code is available in this other repo [RAMON](https://github.com/tiaqui/ramon) (in Spanish), in which you can find functionality for:

* Obtaining input audio in pre-defined cycles of fixed length.
* Applying the microphone's inverse filter.
* Octave band and fractional octave-band filtering.
* Frequency weighting according to the _A_, _C_ or _Z_ curves.
* Linear integration or exponential (_Fast_ or _Slow_) time weightings.
* Recording a calibration tone and saving its RMS level.
* RMS level calculation of the input audio.
* dB<sub>SPL</sub> and L<sub>eq</sub> levels measurement.
* Correction of the octave-band or fractional octave-band levels.
* _DataFrame_, _.npy_ or _HDF5_ format storage of values.
* GUI execution for: 
    - the transformation of values between compressed and tabular formats.
    - _datetime_ index generation.
    - customization of the measured period to preview.
    - modify the temporal integration increasing granularity.
    - percentile calculation.
    - L<sub>den</sub> calculation.
    - visualization of the desired values.
    - storing the information in a tabular file (_.xlsx_ for example).

# Notebooks
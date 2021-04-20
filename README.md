# My engineering thesis

In this repo you'll find Jupyter notebook examples of different funcionality that was used for my **engineering thesis**: **_"Development of a low-cost urban acoustic monitoring station"_**.

## Overview
This thesis was developed to finally obtain my engineering degree at the [UNTREF](https://www.untref.edu.ar/carrera/ingenieria-de-sonido) university. You can see the full _.pdf_ thesis document (in Spanish) [here](doc/Iaquinta,%20Tomás%20-%20Desarrollo%20de%20una%20estación%20de%20monitoreo%20acústico%20urbano%20de%20bajo%20costo.pdf).

What inspired this work was one of the main problems that many experience nowadays (mainly in big cities like BA, where I live): **noise**. The first step in taking action to reduce noise is to have precise objective measurements of it. Traditionally, this was done by professionals using sound level meters. In the present, there are devices (acoustic monitoring stations) that overcome the limitations of that process (mainly the time and space limitations in the measurements).

This thesis objective is to **develop a low-cost device that allows acoustic noise measurements continuously, remotely and autonomously**.

## Developed code
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

## Notebooks

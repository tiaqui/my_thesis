# My engineering thesis

In this repo you'll find Jupyter notebook examples of different funcionality that was used for my **engineering thesis**: **_"Development of a low-cost urban acoustic monitoring station"_**.

## Overview
This thesis was developed to finally obtain my engineering degree at the [UNTREF](https://www.untref.edu.ar/carrera/ingenieria-de-sonido) university. You can see the full _.pdf_ thesis document (in Spanish) [here](doc/Iaquinta,%20Tomás%20-%20Desarrollo%20de%20una%20estación%20de%20monitoreo%20acústico%20urbano%20de%20bajo%20costo.pdf).

What inspired this work was one of the main problems that many experience nowadays (mainly in big cities like BA, where I live): **noise**. The first step in taking action to reduce noise is to have precise objective measurements of it. Traditionally, this was done by professionals using sound level meters. In the present, there are devices (acoustic monitoring stations) that overcome the limitations of that process (mainly the time and space limitations in the measurements).

This thesis objective was to **develop a low-cost device that allows acoustic noise measurements continuously, remotely and autonomously**. Since most of the similar alternatives in the market are cost-prohibitive, since they are usually manufactured with a portable and standard-compliant professional sound level meters.

## Development
The designed device uses a [Raspberry Pi 2 Model B](https://www.raspberrypi.org/products/raspberry-pi-2-model-b/) along with a [digital I2S MEMS microphone](https://www.adafruit.com/product/3421) to capture, process and measure the input sound. A customized housing for the microphone was built using a 3D printer, looking like this:
![3D device model](/img/model.jpg)

For the audio processing the below block diagram was proposed, in which the focus was to adapt the regular processes that traditional analog sound level meters perform to a digital device. ![](_missing image_)

In the notebooks included in this repo we will go over most of the functions and processes neccessary to have a functioning monitoring station.

##  Digital processing - Code
The core code is available in dthis other repo [RAMON](https://github.com/tiaqui/ramon) (in Spanish), in which you can find functionality for:

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

The basic monitoring station's operation can be seen in the image below:

![Monitoring station block diagram](/img/block_diagram.png)

## Notebooks
1. [**Inverse filtering design and implementation**.](/notebooks/1.%20Inverse%20filter%20design%20and%20implementation.ipynb)

![Notebook 1](/img/notebook_1.png)

2. [**A, C and Z frequency weightings**.](/notebooks/2.%20A,%20C%20and%20Z%20frequency%20weightings.ipynb)

![Notebook 2](/img/notebook_2.png)

3. [**Octave and one-third octave band filters**.](/notebooks/3.%20Octave%20and%20one-third%20octave%20band%20filters.ipynb)

![Notebook 3](/img/notebook_3.png)

4. [**Slow and Fast time weightings**.](/notebooks/4.%20Slow%20and%20Fast%20time%20weightings.ipynb)

![Notebook 4](/img/notebook_4.png)

5. [**Acoustic measurements**](/notebooks/5.%20Acoustic%20measurements.ipynb)

![Notebook 5](/img/notebook_5.png)

## Requirements
To run the notebooks the folowing packages are required:
* `numpy`
* `scipy`
* `matplotlib`
* `datetime`
* `pandas`
* `seaborn`
* `IPython`
* `ipywidgets`

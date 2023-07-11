# Implementation of a broad-band GSC beamformer structure

This project consists of the implementation of a generalized sidelobe canceller (GSC) beamformer structure. The GSC structure was used both to steer the array and to provide minimum variance beamforming. This structure represents a time-domain beamformer implementation. The task is to enhance the signal of interest coming from some known directions and suppress the interfering signals at the same time.

This array/beamformer structure is characterized as a multichannel spatial/spectral filter. The objective in broad-band beamformer design is to find a set of weighting coefficients ***w*** , which provide an output having “better” signal-to-noise (SNR) and signal-to-interference (SIR) characteristics than would be observed at the output of a single sensor.

## Description

The [pyroomacoustics library](#references) was used to simulate a room with certain dimensions, under certain reverberation and SNR conditions; a microphone array (uniform linear or uniform circular) and the voice sources located at a certain distance and direction from the array. With this library as core, the [`dataSetAudio.py`](./docs/dataSetAudio.py) module was created, which allows to consider and apply certain functionality and configurations to the simulated acoustic scenario.
As part of this module, a method of the class that performs VAD on the audio files that are passed for the generation of the samples was implemented. [The WebRTCvad Python API](#references) was selected for this task, that was developed by Google for the Webrtc project which aids developers create real time web communication services.

Once the samples at the output of the array are available, the next step is to use the [`matrestriction.py`](./docs/matrestriction.py) module to generate the constraint matrices that constitute the representation space of the voice signals, determine the beamformer coefficients and the blocking matrix of the GSC. Available these, algorithm is carried out.


   ### Dependencies

   - [pyroomacoustics](https://pypi.org/project/pyroomacoustics/) (@LCAV)
   - [webrtcvad](https://pypi.org/project/webrtcvad/) (@wiseman)



## How to use it

1. Import the modules:
   ```py
   import storedata as sd
   import dataSetAudio as dsa
   import matrestriction as mr
   ```

2. Set parameters and retrieve signal modeling matrices.:
   ```python
   #Set parameters for Beamforming
   M = 8         #number of mics in the array
   J = 48        #48 tap delay lines for each channel
   NumAng = 181  #considering 181 DOAs in the scan of directions
       
   # Get the restriction matrices for:
   # -a Beamformer of 8 channels
   # -48 tap delay lines for each channel
   # -scan of 181 DOAs
   # -fs=16000 Hz
   sd.load("restriction_matrices")
   C_list,fr_list,Bk_list = sd.C_list,sd.fr_list,sd.Bk_list # restriction matrices for ULA
   Cc_list,frc_list,Bkc_list = sd.Cc_list,sd.frc_list,sd.Bkc_list # restriction matrices for UCA
   ```

3. Simulate the soundstage:
   ```py
   #Get the audio files and sampling frequency
   fs,data1 = wavfile.read("../voice_data/interfering_voice.wav")
   _,data2 = wavfile.read("../voice_data/desired_voice.wav")

   #Create an object
   direcc_arrival = [24, -41]
   reverb = 0.22
   snr = 0
   distances = [1.7, 2.1]

   ex = dsa.dataSetAudio(direcc_arrival, reverb, snr, distances)
   M_Sigs = ex.SimuData(_, [data1,data2], fs, M, delay=[0,0],N=2,array_mic="ULA",dataset=False)
   ```
4. GSC Algorithm

### Examples

A [notebooks](./docs/GSC.ipynb) with the implementation is available.

A couple of audio files and a binary file with the matrices needed for the case of a beamformer structure consisting of 8 microphones and  48 Tap Delay Lines per channel can be found in the `docs` and `voice_data` folders of the GitHub repository.

## References

[[1]](https://ieeexplore.ieee.org/document/8461310)
   Scheibler, R., Bezzam, E. and Dokmanić, I.: ``Pyroomacoustics: A Python Package for Audio Room Simulation and Array Processing Algorithms'', in: 2018 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), pp. 351–355, doi: 10.1109/ICASSP.2018.8461310.

[[2]](https://pypi.org/project/webrtcvad/)
   Python interface to the Google WebRTC Voice Activity Detector (VAD), Project description.

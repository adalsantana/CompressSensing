# CompressSensing
Github repository for MATLAB code produced for combining compress sensing with the ultrasound process

The goal of this code was to be able to produce a 3D ultrasound image while using only a single sensor. Our physical design consisted of a single transistor with a plastic mask over it with varying height levels and this mask with the induced timing delays and the reflected signal would serve as the sensor array. 
Compress sensing is the mathematics utilized to try recreating a successful image with only the most relevant received timing delays to reconstruct the signal maintaining their weight
This project in its current state did not successfully produce a 3D image

Quick side note: When running the code if you get an error that reads you cannot
set the read-only preoprty 'number_elements' of kWaveTransducer. As of this writing I believe that issue just has to do with the script having been previously ran and you now have a kWaveTransducer object that already exists. To fix it either comment out the part where you create the kwave transducer or just type clear in the command window to erase any transducer elements you may have that already exist

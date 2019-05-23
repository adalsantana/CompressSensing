# CompressSensing
Github repository for MATLAB code produced for combining compress sensing with the ultrasound process

Quick side note: When running the code if you get an error that reads you cannot
set the read-only preoprty 'number_elements' of kWaveTransducer. As of this writing I believe that issue just has to do with the script having been previously ran and you now have a kWaveTransducer object that already exists. To fix it either comment out the part where you create the kwave transducer or just type clear in the command window to erase any transducer elements you may have that already exist

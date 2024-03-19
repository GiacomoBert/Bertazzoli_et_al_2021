To use the scripts for TMS-EEG artifact rejection, please follow the steps below:
1. Make sure EEGLAB and the wavelet toolbox are added to the MATLAB path (the wavelet toolbox is also included in the script folder in case it is not part of your MATLAB package).
2. Add the script folder to the MATLAB path.
3. Prepare the raw continuous EEG data in EEGLAB data structure (for one condition only).
4. Set in cfg.EventCode the event code for the TMS pulses (needs to be a string or a numeric value).
5. Set in cfg.TMSEEGrootFolder the path of the folder to store results.  
6. Use the main function ARTISTMain to preprocess your TMS-EEG data.

An example dataset (exampleDataset.mat) is provided in the script folder for testing purpose. 

Please report bugs by emailing to wwumed@stanford.edu.

Please cite the following paper if the algorithm is used in your publication:
Wu W, Keller C, Rogasch N, Longwell P, Shpigel E, Rolle C, Etkin A. ARTIST: a fully automated artifact rejection algorithm for single-pulse TMS-EEG data. Human Brain Mapping, 2018, DOI: 10.1002/hbm.23938.


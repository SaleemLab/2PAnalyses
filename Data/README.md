# Data
For pre-processing raw data (two-photon, ephys, bonsai, ...) and loading/saving data. To start with, open the file main_preprocess.py and follow the instructions in the comments.

For pre-processing of bouton data (using a different method for registration use main_zregister.npy and follow the instructions there. This will z-register and detect ROIs. After it run use suite2p to curate the folder 'plane'.
After curation use main_preprocess.npy to finalise processing. Do not use a z-stack in this case.


*** This code is cloned from Schroeder Lab (https://github.com/Schroeder-Lab/Data) 
The input and output directories (ritd-ag-project-rd01ie-asale69) have been modified for Saleem Lab.  


Setup for z-motion correction (boutons)
1. Clone/Download this repository 
2. Import conda environment from yaml file
    
    conda env create -f ""C:\PATH-TO-DATA-REPO\Data\refactor.yaml"
    
3. Launch anaconda prompt and activate conda environment 
    
    conda activate refactor 

4. Install Spyder IDE (Python 3.8) within the conda environment 
    
    pip install spyder 
    
5. After spyder has been installed, lanch it in the conda environmet:

    spyder 
    
6. Edit the directories_to_register() function in user_defs.py 
    with the mouse ID, session number and stimulus name 
    
7. You may need to tweek the create_ops_boutton_registration() in user_defs.py 
   to change the number of planes, frame rate, etc. The current parameters 
   are set for boutons with 8 planes, 256x256 pix, 7.28hz and 2 channels. 
   
8. To run the z-motion correction pipeline, run the script in main_zregister.py

 
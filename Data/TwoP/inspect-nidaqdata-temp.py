# -*- coding: utf-8 -*-
"""
Temp functions to read and inspect nidaq files saved by Schnoider lab 
Created on Sat Feb  8 14:37:30 2025

@author: sonali.sriranga
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import random
import sklearn
import seaborn as sns
import scipy as sp
from matplotlib import rc
import matplotlib.ticker as mtick
import matplotlib as mpl
import pandas as pd
import os
import glob
import pickle
import traceback

from Data.TwoP.general import get_ops_file
from Data.TwoP.runners import *
from Data.Bonsai.extract_data import *
from Data.user_defs import *
from Data.Bonsai.behaviour_protocol_functions import stimulus_prcoessing_dictionary

# %%
dirs = define_directories()

csvDir = dirs["dataDefFile"]
s2pDir = dirs["preprocessedDataDir"]
zstackDir = dirs["zstackDir"]
metadataDir = dirs["metadataDir"]


# %% Inspect nidaq bin file
for i in range(len(database)):
    if database.loc[i]["Process"]:
        try:
            print("reading directories" + str(database.loc[i]))
            
            # Get directories for suite2p, zstack, metadata, and save paths
            s2pDirectory, zstackPath, metadataDirectory, saveDirectory = read_csv_produce_directories(
                database.loc[i], s2pDir, zstackDir, metadataDir
            )

            
            # Define the path for the Nidaqinput1.bin file 
            piezoDir = ops["data_path"][0]
            
            planePiezo = get_piezo_data(ops)
   
            # This function reads the nidaq.bin file and the correspoinding channel 
            nidaq, channels, nt = get_nidaq_channels(piezoDir, plot=True)

            nidaq_df = pd.DataFrame(nidaq, columns=channels)
            csv_path = os.path.join(piezoDir, "nidaq-to-inspect.csv")
            nidaq_df.to_csv(csv_path, index=False)
            
        except Exception as e:
            print(f"Error processing session {i}: {e}")




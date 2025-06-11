# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 16:20:30 2025

@author: sonali.sriranga
"""

import sys
sys.path.append("C:/Users/sonali.sriranga/Documents/GitHub")

from suite2p.registration.zalign import compute_zpos
from joblib import Parallel, delayed
import numpy as np
import time
import traceback
import io
import os
import cv2
import skimage.io
import glob
import pickle
import scipy as sp
import tifffile
import re
import warnings
from Data.TwoP.process_tiff import *
from Data.TwoP.preprocess_traces import *
from Data.Bonsai.extract_data import *
from Data.Bonsai.behaviour_protocol_functions import *
from Data.TwoP.general import *
from Data.user_defs import create_2p_processing_ops, directories_to_register
import matplotlib.gridspec as gridspec


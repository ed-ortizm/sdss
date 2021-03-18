#! /usr/bin/env python3
import glob
import os
import time

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from constants_sdss import science_arxive_server_path, spectra_path
from proc_sdss_lib import DownloadData

################################################################################
ti = time.time()
################################################################################
# Data processing

################################################################################

tf = time.time()

print(f'Running time: {tf-ti:.2f} [seg]')

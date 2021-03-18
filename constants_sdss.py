import numpy as np

n_wave = 4_001
wave_master = np.linspace(3_500, 7_500, n_wave)

working_dir = '/home/edgar/oso/spectra_data_curation'
spectra_path = '/home/edgar/oso/sdss_data/spectra'
raw_spectra_data_path = f'{spectra_path}/raw_spectra'
science_arxive_server_path = f'{spectra_path}/sas'
processed_spectra_path = f'{spectra_path}/processed_spectra'

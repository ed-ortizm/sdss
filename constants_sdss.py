import numpy as np

n_wave = 4_001
wave_master = np.linspace(3_500, 7_500, n_wl)

spectra_data_path = '/home/edgar/zorro/sdss_data/spectra_sdss'
raw_spectra_data_path = f'{spectra_data_path}/raw_spectra'
science_arxive_server_path = f'{spectra_data_path}/sas'
processed_spectra_path = f'{spectra_data_path}/processed_spectra'

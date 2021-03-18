import numpy as np

n_waves = 4_001
wave_master = np.linspace(3_500, 7_500, n_waves)

working_dir = '/home/edgar/oso/spectra_data_curation'
spectra_path = '/home/edgar/oso/sdss_data/spectra'
science_arxive_server_path = f'{spectra_path}/sas'
processed_spectra_path = f'{spectra_path}/processed_spectra'

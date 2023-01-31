<!-- # SDSSdata
SDSS spectroscopy data
Project Name & Description: Give your project a name and a brief, clear, and concise description of what the library does.

Installation: Explain how to install the library, including any necessary dependencies.

Usage: Provide examples of how to use the library, including code snippets.

API Reference: If your library has a public API, document it in detail, including function signatures and descriptions.

Contributing: Provide guidelines for contributing to the project, including how to report bugs and request features.

License: Include information on the license under which the project is released.

Acknowledgements: Give credit to any third-party libraries used in the project.

Further Reading: Provide links to additional resources, such as tutorials or documentation, for those who want to learn more. -->

# SDSS Spectra Library

A library for downloading and processing raw spectra from the Sloan Digital Sky Survey (SDSS) DR16.
This library leverages the power of the multiprocessing library from python to run most of the data processing tasks using CPU parallelization.

## Installation

1. Clone the repository
```bash
git clone https://github.com/ed-ortizm/sdss
```
2. Install dependencies

```python
pip install -r requirements.txt
```

3. Install the library
```python
python setup.py install
```

## Modules & library structure

```csharp
src
└── sdss
    ├── describe.py
    ├── download.py
    ├── metadata.py
    ├── process
    │   ├── deredspectra.py
    │   ├── filter.py
    │   ├── indefinite_values.py
    │   ├── inputting.py
    │   ├── interpolate.py
    │   └── sample.py
    ├── raw
    │   ├── data.py
    └── utils
        ├── configfile.py
        ├── managefiles.py
        ├── parallel.py
        └── timer.py
```
# Usage
For examples of how to use of the scripts found in,
```csharp
.
├── download.ini
├── download.py
├── process
│   ├── inputting.ini
│   ├── inputting.py
│   ├── interpolate.ini
│   ├── interpolate.py
│   ├── remove.ini
│   └── remove.py
├── raw
│   ├── ebv_values.ini
│   ├── ebv_values.py
│   ├── raw.ini
│   └── raw.py
└── sample
    ├── describe.ini
    ├── describe.py
    ├── explore_bins.ipynb
    ├── sample.ini
    ├── sample.py
    ├── train_sets.ini
    └── train_sets.py
```

please notice that each has a companion configuration file (e.g. download.ini, raw.ini, process.ini) where you can specify the parameters.
For a description of each script please check the README.md file of each directory.

# License

This library is licensed under the MIT License. The MIT License is a permissive free software license that allows for reuse within proprietary software provided all copies of the licensed software include a copy of the MIT License terms and the copyright notice.

The full text of the license can be found in the LICENSE file.

# Acknowledgements
This library was built with the following Python packages:

* Astropy
* Matplotlib
* NumPy
* pandas
* SciPy
* setuptools
* sfdmap

I thank the developers of these packages for their efforts in making the scientific Python ecosystem so powerful.
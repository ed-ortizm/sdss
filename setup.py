from setuptools import setup, find_packages

with open("README.md", "r") as file:
    long_description = file.read()

setup(
    name="sdss",
    version="0.0.1",
    author="Edgar Ortiz",
    author_email="ed.ortizm@gmail.com",
    packages = find_packages(
        where='src',
        include=['[a-z]*'],
        exclude=['old_code']
    ),
    package_dir = {"":"src"},
    # packages=["autoencoders"],
    # package_dir = {"autoencoders":"src/autoencoders"},
    description="Module to process sdss spectra for outlier detection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ed-ortizm/spectra-processing",
    license="MIT",
    keywords="astrophysics, galaxy, Machine Learning, SDSS"
)

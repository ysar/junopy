
import setuptools
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="JunoPy",
    version="0.0.1",
    author="Yash Sarkango",
    author_email="sarkango@umich.edu",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    description=
        "A python package for analysis of data from the Juno Spacecraft.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ysar/JunoPy",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6'
)

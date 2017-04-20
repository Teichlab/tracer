import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name='tracer',
    version=0.1,
    author="Mike Stubbington",
    entry_points={
        'console_scripts': [
            'tracer = tracerlib.launcher:launch'
        ]
    },
    author_email="mike.stubbington@sanger.ac.uk",
    description="Reconstruction of T-Cell receptor sequences from single-cell RNA-seq data",
    licence="Apache",
    keywords="biopython genetics",
    url="https://github.com/teichlab/tracer",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.66",
        "cycler>=0.10.0",
        "decorator>=4.0.9",
        "matplotlib>=1.5.1",
        "networkx>=1.11",
        "numpy>=1.11.0",
        "pandas>=0.18.0",
        "prettytable>=0.7.2",
        "pydotplus>=2.0.2",
        "pyparsing>=2.0.3",
        "python-dateutil>=2.5.2",
        "python-Levenshtein>=0.12.0",
        "pytz>=2016.3",
        "scipy>=0.17.0",
        "seaborn>=0.7.0",
        "six>=1.10.0",
        "mock>=2.0.0",
        "future>=0.15.2"
    ]
)

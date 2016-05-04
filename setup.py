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
            'tracer=tracer.tracerlib.__main__'
        ]
    },
    author_email="mike.stubbington@sanger.ac.uk",
    description="Reconstruction of T-Cell receptor sequences from single-cell RNA-seq data",
    licence="Apache",
    keywords="biopython genetics",
    url="https://github.com/teichlab/tracer",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "cycler",
        "decorator",
        "matplotlib",
        "networkx",
        "numpy",
        "pandas",
        "prettytable",
        "pydotplus",
        "pyparsing",
        "python-dateutil",
        "python-Levenshtein",
        "pytz",
        "scipy",
        "seaborn",
        "six",
        "future"
    ]
)
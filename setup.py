import os
import glob
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

#extra_files = package_files('test_data')

def get_requirements():
    with open("/requirements.txt", "rt", encoding="utf-8") as fh:
        return [line.strip() for line in fh.readlines()]

setup(
    name='tracer',
    version=0.6,
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
    #package_data={'tracer': extra_files},
    install_requires=get_requirements()
)

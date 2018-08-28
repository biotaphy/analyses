# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='ancestral_reconstruction',
    version='0.1.0',
    description='Biotaphy package for generating ancestral reconstructions',
    long_description=readme,
    author='Biotaphy Team',
    author_email='cjgrady@ku.edu',
    url='https://github.com/biotaphy/ancestral_state',
    license=license,
    packages=find_packages(exclude=('tests', 'docs', 'sample_data')),
    scripts=['bin/ancestral_distribution.py', 'bin/ancestral_state.py'],
    install_requires=[
        'dendropy>=4.0.0', 
        'matplotlib', 
        'numpy>=1.11.0', 
        'scipy>=1.0.0']
)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from setuptools import setup, find_packages

with open("README.md") as f:
    long_description = f.read()

with open('pygalaxy/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__author__'):
            author = line.split('=')[-1]
        if line.startswith('__version__'):
            version = line.split('=')[-1]

setup(
    name='pygalaxy',
    version=version,
    description='',
    long_description=long_description,
    author=author,
    url='https://github.com/correac/pygalaxy',
    license="BSD",
    keywords=['pygalaxy', 'cosmology', 'NFW', 'hot halo', 'accretion'],
    classifiers=['Development Status :: 0',
                 'Intended Audience :: Developers',
                 'Natural Language :: English',
                 'Programming Language :: Python :: 3.0'],
    install_requires=['numpy',
                      'scipy',
                      'h5py',
                      'matplotlib']
      #entry_points={
      #'console_scripts': [
      #      'mahstery = pygalaxy.pygalaxy:run'
      #  ],
      #},
      #packages=find_packages()
)

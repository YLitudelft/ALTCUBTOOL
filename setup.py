# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 10:39:46 2020

@author: lyh
"""

from setuptools import setup, find_packages

setup(
    name='ALTCUBTOOL',
    packages=find_packages(),
    version='0.0.0',
    install_requires=[         
        'matplotlib', 'ocs_io', 'pickle', 'os'
    ]
)
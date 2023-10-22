#!/usr/bin/env python3
from __future__ import print_function, absolute_import
from setuptools import setup, find_packages
from distutils.sysconfig import get_python_lib
import os
import glob
import sys
pwd_ = os.path.dirname(os.path.abspath(__file__))

print("===================== IGmodel Setup.py =========================")
print("The setup.py script should be executed from the build directory.")
print("pwd: " + pwd_)

setup(
    name='IGmodel',
    version='0.0.1',
    packages=find_packages(),
    description='General property calculator model for an ideal gas mixture',
    author='Micah Best',
    platforms=["Linux"],
    build_requires=['numpy>=1.8.0', 'setuptools', 'pylint'],
    install_requires=['numpy>=1.8.0', 'scipy>=0.12.0', 'pylint'],
    license='MIT',
    author_email='micahbest21@gmail.com',
    keywords='ideal gas model',
)
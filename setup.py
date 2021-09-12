# -*- encoding: utf-8 -*-
'''
Description:  PyPI     
@created   : 2021/09/08 09:33:59
'''

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="raser",
    version="2.1.1",
    author="Xin Shi",
    author_email="xin.shi@outlook.com",
    description="SiC Timing Detector Simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dt-np/raser",
    packages=find_packages(),
    classifiers=[
				"Programming Language :: Python :: 3",
				"License :: OSI Approved :: MIT License",
				"Operating System :: OS Independent",
    			]
	)
